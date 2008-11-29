"""Read a file with sequence features and create a file with sparse real features based on word distance histogram."""

from MGT.Shogun.Util import *
from shogun.Features import *
from shogun.Classifier import *
from shogun.Distance import *
from MGT.PredProcessor import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-feat",
        action="append", type="string",dest="inFeat"),
        make_option("-c", "--svm-penalty",
        action="store", type="float",dest="C",default=1.),
        make_option("-m", "--mode",
        action="store", type="choice",choices=("train","test","predict"),dest="mode",default="train"),
        make_option("-e", "--method",
        action="store", type="choice",choices=("svm","knn","knn-svm"),dest="method",default="svm"),
        make_option("-k", "--knn-k",
        action="store", type="int",dest="knnK",default=5),
        make_option("-l", "--train-label",
        action="store", type="string",dest="trainLabel",default="-1"),
        make_option("-t", "--predict-thresh",
        action="store", type="string",dest="predThresh",default="-2000"),
        make_option("-o", "--model-name",
        action="store", type="string",dest="modelRoot",default="mod"),
        make_option("-p", "--pred-file",
        action="store", type="string",dest="predFile"),
        make_option("-a", "--lab-to-id",
        action="store", type="string",dest="labToId"),
        make_option("-g", "--gos-dir",
        action="store", type="string",dest="gosDir",default=None),
        make_option("-n", "--check-null-hyp",
        action="store_true", dest="checkNull",default=False),
        make_option("-s", "--svm-srm",
        action="store_true", dest="useSrm",default=False),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


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



def analyzePerf(labTest,labPred,labToId,name):
    pm = perfMetrics(labTest,labPred,balanceCounts=False)
    pm.setLabToId(labToId)
    print pm.toIdSpeStr()
    print pm.toIdSenStr()
    pm.confMatrCsv(name+'cm.csv')


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
        assert not labSeen[0]
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
    svmMul = SvmOneVsAll(maxLabel=maxLab)
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

class App:

    def __init__(self,opt):
        self.opt = opt
        feat = []
        lab = []
        inFeat = [ l for l in opt.inFeat ]
        if len(inFeat) < 3:
            inLast = inFeat[-1]
            for i in range(len(inFeat),3):
                inFeat.append(inLast)
        print inFeat
        self.inFeat = inFeat
        for i in xrange(len(inFeat)):
            # just append the previous loaded object if the file name is the same
            if not (i > 0 and inFeat[i] == inFeat[i-1]):
                f = SparseRealFeatures()
                l = f.load_svmlight_file(inFeat[i]).get_labels().astype('i4')
            feat.append(f)
            lab.append(l)
        self.feat = feat
        self.lab = lab

    def processOptions(self):
        feat = self.feat
        lab = self.lab
        opt = self.opt
        
        trainLabels = [ int(l) for l in opt.trainLabel.split(',') ]

        thresh = [ float(t) for t in opt.predThresh.split(',') ]
        if len(thresh) == 3:
            thresh = numpy.arange(thresh[0],thresh[1],thresh[2])
        elif len(thresh) == 1:
            thresh = numpy.array(thresh,dtype='f4')
        else:
            raise AssertionError("threshold must be 'start,end,step' or 'value'")


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
            if hasattr(opt,"labToId") and opt.labToId is not None:
                labToId = loadObj(opt.labToId)
            else:
                labToId = n.arange(20000)
            self.labToId = labToId
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
                for t in thresh:
                    labPred = svmMul.classify(opt=Struct(thresh=t,useSrm=opt.useSrm))
                    print "Threshold %.3f" % t
                    yield labPred
            elif opt.method == "knn":
                dist = SparseEuclidianDistance(feat[0],feat[1])
                mod = KNN(opt.knnK,dist,Labels(lab[0].astype('f8')))
                mod.train()
                labPred = mod.classify().get_labels()
                yield labPred
            elif opt.method == "knn-svm":
                dist = SparseEuclidianDistance(feat[0],feat[2])
                knn = KNN(opt.knnK,dist,Labels(lab[0].astype('f8')))
                knn.train()
                n_test = feat[2].get_num_vectors()
                ind_neighb = numpy.zeros((n_test,opt.knnK),dtype='i4')
                print "Computing KNN list..."
                knn.get_neighbours(ind_neighb)
                labPred = numpy.zeros(n_test,dtype='i4')
                print "Training neighbours' SVMs..."
                for iTest in xrange(n_test):
                    svmTrFeat = feat[0].subsample(ind_neighb[iTest])
                    svmTrLab = lab[0][ind_neighb[iTest]]
                    if (svmTrLab == svmTrLab[0]).all():
                        labPred[iTest] = svmTrLab[0]
                        if iTest % 100 == 0:
                            print "All k neighbours have one label %i for samp %i" % (labPred[iTest],iTest)
                    else:
                        svmTsFeat = feat[2].subsample(numpy.asarray([iTest],dtype='i4'))
                        labPred[iTest] = svmOneVsAllOneStep(feat=(svmTrFeat,svmTsFeat),
                                lab=(svmTrLab,),
                                opt=Struct(C=opt.C,thresh=thresh[0],useSrm=False))
                        if iTest % 100 == 0:
                            print "SVM selected label %i from %s for samp %i" % (labPred[iTest],svmTrLab,iTest)

                yield labPred

    def run(self):
        opt = self.opt
        for labPred in self.processOptions():
            if opt.mode == "test":
                analyzePerf(self.lab[2],labPred,self.labToId,"mod")
            elif opt.mode == "predict":
                print self.labToId
                print numpy.bincount(labPred.astype('i4'))
                dumpObj(labPred,opt.predFile)
                if opt.gosDir is not None:
                    analyzePredGos(labPred,opt.gosDir)


opt,args = getProgOptions()

app = App(opt)

app.run()

