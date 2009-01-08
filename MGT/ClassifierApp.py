"""Application-level interface to supervised classifiers."""

from MGT.Shogun import *
from MGT.Svm import *
from shogun.Features import *
from shogun.Classifier import *
from shogun.Distance import *
from MGT.PredProcessor import *
from MGT.App import *


def softBinDec(y):
    return numpy.sign(y)*(1-numpy.exp(-numpy.abs(y)))

def srm(wf,c,l,y):
    a = 1-l*y
    return numpy.exp(-(wf + c*numpy.select([a>0],[a],default=0).sum())/(c*len(y)))

def srmBinDec(wf,c,l,y):
    rm = srm(wf,c,l,y)
    sy = softBinDec(y)
    return 

class ModelLinIO:
    
    def getData(self,svm):
        d = Struct()
        d.bias = svm.get_bias()
        d.bias_enabled = svm.get_bias_enabled()
        d.w = svm.get_w()
        return d

    def save(self,out,svm):
        d = self.getData(svm)
        dumpObj(d,out)

    def load(self,inp,svm):
        d = loadObj(inp)
        svm.set_bias(d.bias)
        svm.set_bias_enabled(d.bias_enabled)
        svm.set_w(d.w)

def svmOcasFact():
    return SVMOcas(SVM_OCAS)


class ModFact:

    def loadModels(self,labels=None):
        if labels == None:
            labels = self.getLabels()
        for label in labels:
            yield (label,self.loadMod(label))

    def getMaxLabel(self):
        labels = self.getLabels()
        return labels[-1]

    def fromOther(self,other,labels=None):
        for label,mod in other.loadModels(labels=labels):
            self.saveMod(label,mod)

class ModFileFact(ModFact):

    maxLabel = 100000
    def __init__(self,modelRoot,svmFact):
        assert '-' not in modelRoot,"model name cannot contain dash symbol"
        self.modelRoot=modelRoot
        self.numLabDig = int(numpy.ceil(numpy.log10(self.maxLabel))+1)
        self.svmFact = svmFact

    def name(self,label):
        labFormat = "%0"+str(self.numLabDig)+"i"
        return os.path.join(self.modelRoot,labFormat % label)

    def saveMod(self,label,svm):
        fileName = self.name(label)
        makeFilePath(fileName)
        ModelLinIO().save(fileName,svm)

    def loadMod(self,label):
        svm = self.svmFact()
        ModelLinIO().load(self.name(label),svm)
        return svm

    def getLabels(self):
        return sorted([ int(modFile) for modFile in os.listdir(self.modelRoot) ])



class ModMemFact(ModFact):

    def __init__(self):
        self.svms = {}

    def saveMod(self,label,svm):
        self.svms[label] = svm

    def loadMod(self,label):
        return self.svms[label]

    def getLabels(self):
        return sorted(self.svms.keys())





class SvmOneVsAll:

    def __init__(self,maxLabel=None):
        self.maxLabel = maxLabel
        #keep bin svms in memory by default
        self.modFact = ModMemFact() 

    def setFeat(self,feat):
        self.feat = feat

    def setLab(self,lab):
        self.lab = lab
        if self.maxLabel is None and self.lab is not None:
            self.maxLabel = max(self.lab)

    def setBinSvm(self,modFact):
        self.modFact = modFact

    def classifyBin(self):
        modIter = self.modFact.loadModels()
        feat = self.feat
        nSamp = feat.get_num_vectors()
        binLab = numpy.zeros((nSamp,self.maxLabel+1),dtype='f8')
        labSeen = numpy.zeros(binLab.shape[1],dtype=bool)
        labShog = Labels(numpy.zeros(nSamp,dtype='f8'))
        for (label,svm) in modIter:
            svm.set_features(feat)
            w = svm.get_w()
            n_w = len(w)
            if feat.get_num_features() < n_w:
                feat.set_num_features(n_w)
            labPred = svm.classify(labShog).get_labels()
            binLab[:,label] = labPred
            labSeen[label] = True
        #assert not labSeen[0] #why did I need this?
        self.binLab = binLab
        self.labSeen = labSeen
        print "SvmOneVsAll.classifyBin() did not see models for labels: "+','.join(["%s" % l for l in n.where(labSeen == False)[0]])

    def trainBin(self,trainLabel,opt):
        feat = self.feat
        lab = self.lab
        labBin = numpy.select([lab == trainLabel], [1.], default=-1.)
        epsilon=1e-2
        num_threads=1
        svm=SVMOcas(opt.C, feat, Labels(labBin)) #SVMOcas LibLinear (broken, inverts the sign for some datasets)
        #svm.set_C(opt.C*float((labBin==-1).sum())/(labBin==1).sum(),opt.C)
        svm.set_epsilon(epsilon)
        svm.parallel.set_num_threads(num_threads)
        svm.set_bias_enabled(True)
        print "Training label %i ..." % trainLabel
        svm.train()
        return svm

    def trainMany(self,trainLabels=None,opt=Struct(C=1.)):
        if trainLabels is None:
            trainLabels = numpy.unique(self.lab)
        modFact = self.modFact
        for trainLabel in trainLabels:
            svm = self.trainBin(trainLabel,opt)
            modFact.saveMod(trainLabel,svm)

    def computeSrm(self):
        labels = self.lab
        modIter = self.modFact.loadModels()
        self.srm = numpy.zeros(len(self.labSeen),dtype='f8')
        for (label,svm) in modIter:
            w = svm.get_w()
            wf = numpy.inner(w,w)/2.
            c = svm.get_C1()
            assert self.labSeen[label]
            y = self.binLab[:,label]
            l = numpy.select([labels == label],[1.],default=-1.)
            self.srm[label] = srm(wf,c,l,y)

    def getSrm(self):
        return self.srm

    def setSrm(self,srm):
        self.srm = srm

    def classify(self,opt):
        binLab = self.binLab.copy()
        if opt.useSrm:
            for label in numpy.where(self.labSeen)[0]:
                binLab[:,label] = softBinDec(binLab[:,label]) * self.srm[label]
        binLab[:,-self.labSeen] = -1e10
        labMaxBin = binLab.argmax(1)
        rowInd = numpy.arange(len(binLab))
        pred = numpy.select([binLab[rowInd,labMaxBin] >= opt.thresh], [labMaxBin], default=0)
        return pred

def svmOneVsAllOneStep(feat,lab,opt=Struct(C=1.,thresh=-200,useSrm=False)):
    svmMul = SvmOneVsAll()
    svmMul.setLab(lab[0])
    svmMul.setFeat(feat[0])
    svmMul.trainMany(opt=opt)
    svmMul.setLab(None)
    svmMul.setFeat(feat[1])
    svmMul.classifyBin()
    labPred = svmMul.classify(opt=opt)
    return labPred


def analyzePredGos(pred,dataDir):
    meta = loadObj(os.path.join(dataDir,"inp.meta"))
    samp = meta.samp
    assert len(pred) == len(samp)
    idToInd = {}
    for iSamp in xrange(len(samp)):
        idToInd[samp[iSamp]['id']] = iSamp
    assert not idToInd.has_key(0)
    nMated = 0
    nMatedSame = 0
    for iSamp in xrange(len(samp)):
        if pred[iSamp] != 0:
            id_mate = samp[iSamp]['id_mate']
            if idToInd.has_key(id_mate):
                iSampMate = idToInd[id_mate]
                if iSamp < iSampMate and pred[iSampMate] != 0:
                    nMated += 1
                    if pred[iSamp] == pred[iSampMate]:
                        nMatedSame += 1
    print "nMated = %i  nMatedSame = %i nMatedSame/nMated = %.2f" % \
            (nMated,nMatedSame,nMatedSame*100./nMated)
    counts = numpy.bincount(pred)
    print "Assignment ratios: ", counts*100./counts.sum()
    counts[0] = 0
    print "Assignment ratios excl reject: ", counts*100./counts.sum()




class ClassifierApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        optChoicesMode = ("train","test","predict")
        option_list = [
            make_option("-i", "--in-feat",
            action="append", type="string",dest="inFeat"),
            make_option("-c", "--svm-penalty",
            action="store", type="float",dest="C",default=1.),
            make_option("-m", "--mode",
            action="store", type="choice",choices=optChoicesMode,
            dest="mode",default="train",help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            make_option("-e", "--method",
            action="store", type="choice",choices=("svm","knn","knn-svm"),dest="method",default="svm"),
            make_option("-k", "--knn-k",
            action="store", type="int",dest="knnK",default=5,help="Number of nearest neighbours in KNN algorithm"),
            make_option("-d", "--knn-max-dist",
            action="store", type="float",dest="knnMaxDist",default=None),
            make_option("-l", "--train-label",
            action="store", type="string",dest="trainLabel",default="-1"),
            make_option("-t", "--predict-thresh",
            action="store", type="string",dest="thresh",default="-2000"),
            make_option("-o", "--model-name",
            action="store", type="string",dest="modelRoot",default="mod"),
            make_option("-p", "--pred-file",
            action="store", type="string",dest="predFile",default=None),
            make_option("-r", "--perf-file",
            action="store", type="string",dest="perfFile"),
            make_option("-a", "--labels",
            action="append", type="string",dest="labels"),
            make_option("-u", "--lab-unclass",
            action="store", type="int",dest="labUnclass",default=0),
            make_option("-g", "--gos-dir",
            action="store", type="string",dest="gosDir",default=None),
            make_option("-n", "--check-null-hyp",
            action="store_true", dest="checkNull",default=False),
            make_option("-s", "--svm-srm",
            action="store_true", dest="useSrm",default=False),
            make_option(None, "--export-conf-matr",
            action="store_true", dest="exportConfMatr",default=False),
            make_option(None, "--save-conf-matr",
            action="store_true", dest="saveConfMatr",default=False),
        ]
        return Struct(usage = "Construct several types of classifiers or use previously"+\
                " trained classifier for testing and prediction\n"+\
                "%prog [options]",option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        thresh = [ float(t) for t in options.thresh.split(',') ]
        if len(thresh) == 3:
            thresh = numpy.arange(thresh[0],thresh[1],thresh[2])
        elif len(thresh) == 1:
            thresh = numpy.array(thresh,dtype='f4')
        else:
            parser.error("--predict-thresh value must be 'start,end,step' or 'value', received %s" % options.thresh)
        options.thresh = thresh

    def loadFeatures(self):
        opt = self.opt
        labels = opt.labels
        if labels is None:
            labels = [ lf+".idlab" for lf in opt.inFeat ]
        idLab = loadIdLabelsMany(fileNames=labels)
        self.idLab = idLab
        assert len(opt.inFeat) > 0
        dataFeat = loadSparseSeqsMany(opt.inFeat,idLab=idLab)
        idLabSplits = idLab.getSplits()
        feat = []
        lab = []
        data = []
        for split in sorted(idLabSplits):
            idl = idLabSplits[split]
            splitData = idl.selData(dataFeat)
            splitFeat = convSparseToShog(data=splitData,delFeature=True)
            data.append(splitData)
            feat.append(splitFeat)
            lab.append(splitData["label"].astype('i4'))
        # repeat the last split if less than 3
        if len(feat) < 3:
            iLast = len(feat) - 1
            for i in range(len(feat),3):
                data.append(data[iLast])
                feat.append(feat[iLast])
                lab.append(lab[iLast])
        print opt.inFeat
        self.inFeat = opt.inFeat
        self.feat = feat
        self.lab = lab
        self.data = data


    def processOptions(self):
        feat = self.feat
        lab = self.lab
        opt = self.opt
        
        trainLabels = [ int(l) for l in opt.trainLabel.split(',') ]

        thresh = opt.thresh

        if opt.checkNull:
            if opt.mode == "train":
                nrnd.shuffle(lab[0])
            elif opt.mode == "test":
                nrnd.shuffle(lab[2])
            else:
                raise "Null hypothesis checking is undefined for this mode: ",opt.mode

        maxLab = max([ max(l) for l in lab ])
        modFact = ModFileFact(opt.modelRoot,svmOcasFact)
        if opt.mode == "train":
            trainLabCnt = numpy.bincount(lab[0])
            trainLabelsPos = numpy.where(trainLabCnt>0)[0]
            if len(trainLabels) == 1 and trainLabels[0] == -1:
                trainLabels = trainLabelsPos
            else:
                trainLabels = [ trLab for trLab in trainLabels if trLab in trainLabelsPos ]
            if opt.method == "svm":
                svmMul = SvmOneVsAll(maxLabel=maxLab)
                svmMul.setLab(lab[0])
                svmMul.setFeat(feat[0])
                svmMul.setBinSvm(modFact)
                svmMul.trainMany(trainLabels=trainLabels,opt=opt)

        elif opt.mode in ("test","predict"):
            if opt.method == "svm":
                if len(trainLabels) == 1 and trainLabels[0] == -1:
                    labLoad = None
                    maxLabel = modFact.getMaxLabel()
                else:
                    labLoad = trainLabels
                    maxLabel = max(trainLabels)
                svms = ModMemFact()
                svms.fromOther(modFact,labels=labLoad)
                svmMul = SvmOneVsAll(maxLabel=maxLabel)
                svmMul.setBinSvm(svms)
                if opt.useSrm:
                    svmMul.setFeat(feat[1])
                    svmMul.classifyBin()
                    svmMul.setLab(lab[1])
                    svmMul.computeSrm()
                    srm = svmMul.getSrm()
                    #srm[:] = 1.
                    #svmMul.setSrm(srm)
                    print "SRM = %s" % (svmMul.getSrm(),)
                svmMul.setLab(None)
                svmMul.setFeat(feat[2])
                svmMul.classifyBin()
                labPred = n.zeros((len(thresh),len(lab[2])),dtype='i4')
                for iThresh in xrange(len(thresh)):
                    t = thresh[iThresh]
                    labPred[iThresh] = svmMul.classify(opt=Struct(thresh=t,useSrm=opt.useSrm))
                    print "Threshold %.3f" % t
                return Struct(labPred=labPred,param=n.rec.fromarrays([thresh],names="thresh"))
            elif opt.method == "knn":
                dist = SparseEuclidianDistance(feat[0],feat[1])
                mod = KNN(opt.knnK,dist,Labels(lab[0].astype('f8')))
                if opt.knnMaxDist is not None:
                    mod.set_max_dist(opt.knnMaxDist)
                mod.train()
                labPred = mod.classify().get_labels()
                labUnclass = mod.get_unclass_label()
                labPred[labPred==labUnclass] = opt.labUnclass
                labPred.shape = (1,len(labPred))
                return Struct(labPred=labPred,param=None)
            elif opt.method == "knn-svm":
                assert len(thresh) == 1,"multiple SVM decision thresholds not implemented for knn-svm"
                dist = SparseEuclidianDistance(feat[0],feat[2])
                knn = KNN(opt.knnK,dist,Labels(lab[0].astype('f8')))
                knn.train()
                n_test = feat[2].get_num_vectors()
                ind_neighb = numpy.zeros((n_test,opt.knnK),dtype='i4')
                dist_neighb = numpy.zeros((n_test,opt.knnK),dtype='f8')
                print "Computing KNN list..."
                knn.get_neighbours(ind_neighb,dist_neighb)
                labPred = numpy.zeros(n_test,dtype='i4')
                print "Training neighbours' SVMs..."
                for iTest in xrange(n_test):
                    samp_ind_neighb = ind_neighb[iTest]
                    samp_dist_neighb = dist_neighb[iTest]
                    if opt.knnMaxDist is not None:
                        samp_in_dist = samp_dist_neighb < opt.knnMaxDist
                        samp_ind_neighb = samp_ind_neighb[samp_in_dist]
                        samp_dist_neighb = samp_dist_neighb[samp_in_dist]
                    if len(samp_ind_neighb) > 0:
                        svmTrFeat = feat[0].subsample(samp_ind_neighb)
                        svmTrLab = lab[0][samp_ind_neighb]
                        if (svmTrLab == svmTrLab[0]).all():
                            labPred[iTest] = svmTrLab[0]
                            if iTest % 100 == 0:
                                print "All %s neighbours have one label %i for samp %i" % (len(samp_ind_neighb),labPred[iTest],iTest)
                        else:
                            svmTsFeat = feat[2].subsample(numpy.asarray([iTest],dtype='i4'))
                            labPred[iTest] = svmOneVsAllOneStep(feat=(svmTrFeat,svmTsFeat),
                                    lab=(svmTrLab,),
                                    opt=Struct(C=opt.C,thresh=thresh[0],useSrm=False))
                            if iTest % 100 == 0:
                                print "SVM selected label %i from %s for samp %i" % (labPred[iTest],svmTrLab,iTest)
                    else:
                        labPred[iTest] = opt.labUnclass
                        if iTest % 100 == 0:
                            print "No training samples are within cutoff distance found for samp %i" % (iTest,)

                labPred.shape = (1,len(labPred))
                return Struct(labPred=labPred,param=None)

    def doWork(self,**kw):
        opt = self.opt
        print "Classifier options: \n%s\n" % opt
        assert not (opt.mode == "test" and opt.perfFile is None)
        assert not (opt.mode == "predict" and opt.predFile is None)
        self.loadFeatures()
        pred = self.processOptions()
        if pred is not None:
            idPred = selFieldsArray(self.data[2],["id","label"])
            assert len(idPred) == len(pred.labPred[0])
            pred = Predictions(labPred=pred.labPred,param=pred.param,idPred=idPred)
            if opt.predFile is not None:
                dumpObj(pred,opt.predFile)
            if opt.mode == "test":
                perf = pred.calcPerfMetrics(idLab=self.idLab,
                        confMatrFileStem=opt.perfFile if opt.exportConfMatr else None,
                        keepConfMatr=opt.saveConfMatr)
                dumpObj(perf,opt.perfFile)
            elif opt.mode == "predict":
                print self.labels.labNames()
                print n.bincount(pred.labPred[0].astype('i4'))
                if opt.gosDir is not None:
                    analyzePredGos(labPred,opt.gosDir)



if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ClassifierApp)

