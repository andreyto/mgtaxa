"""Application-level interface to supervised classifiers."""

from MGT.Shogun import *
from MGT.Svm import *
from shogun.Features import *
from shogun.Classifier import *
from shogun.Kernel import SparseGaussianKernel
from shogun.Distance import *
from MGT.PredProcessor import *
from MGT.SvmMulti import *
from MGT.App import *

__all__ = ["ClassifierApp"] 


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
            make_option(None, "--kernel",
            action="store", type="choice",choices=("lin","rbf"),dest="kernel",default="lin"),
            make_option(None, "--rbf-width",
            action="store", type="float",dest="rbfWidth",default=1.),
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
            make_option(None, "--balance-test-counts",
            action="store_true", dest="balanceTestCounts",default=False,
            help="Balance test sample counts across labels when calculating confusion matrix"),
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
        assert not isinstance(labels,str),"Need a sequence of label files"
        if labels is None:
            labels = [ lf+".idlab" for lf in opt.inFeat ]
        idLab = loadIdLabelsMany(fileNames=labels)
        #DEBUG:
        #idLab = idLab.balance(100)
        self.idLab = idLab
        assert not isinstance(opt.inFeat,str),"Need a sequence of feature files"
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
        if "svm" in opt.method:
            if opt.kernel == "lin":
                svmFact = SvmShogLin.factory(C=opt.C)
            elif opt.kernel == "rbf":
                kernel=SparseGaussianKernel(100,opt.rbfWidth)
                #kernel must know its lhs for classification, and
                #it must be the same as it was for training
                kernel.init(feat[0],feat[0])
                svmFact = SvmShogKern.factory(C=opt.C,kernel=kernel)
            if opt.method == "svm":
                modStore = SvmModFileStore(opt.modelRoot,svmFact=svmFact)
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
                svmMul.setSvmStore(modStore)
                svmMul.trainMany(trainLabels=trainLabels)

        elif opt.mode in ("test","predict"):
            if opt.method == "svm":
                if len(trainLabels) == 1 and trainLabels[0] == -1:
                    labLoad = None
                    maxLabel = modStore.getMaxLabel()
                else:
                    labLoad = trainLabels
                    maxLabel = max(trainLabels)
                svms = SvmModMemStore(svmFact)
                svms.fromOther(modStore,labels=labLoad)
                svmMul = SvmOneVsAll(maxLabel=maxLabel)
                svmMul.setSvmStore(svms)
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
                    labPred[iThresh] = svmMul.classify(thresh=t,useSrm=opt.useSrm)
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
                        keepConfMatr=opt.saveConfMatr,balanceCounts=opt.balanceTestCounts)
                dumpObj(perf,opt.perfFile)
                perf.exportMetricsCsv(names=("senMin","speMin","senMean","speMean","acc"),
                        out=stripPathSfx(opt.perfFile)+".csv")
            elif opt.mode == "predict":
                #print n.bincount(pred.labPred[0].astype('i4'))
                #if opt.gosDir is not None:
                #    analyzePredGos(labPred,opt.gosDir)
                pass



if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ClassifierApp)

