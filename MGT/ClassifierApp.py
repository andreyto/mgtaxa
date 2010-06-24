### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


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
    
    batchDepModes = ("trainScatter",)
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        optChoicesMode = ("trainScatter","train","test","predict")
        option_list = [
            make_option("-i", "--in-feat",
            action="append", type="string",dest="inFeat"),
            make_option(None, "--in-feat-format",
            action="store", type="choice",choices=featIOFormats,
            dest="inFeatFormat",default=defFeatIOFormat),
            make_option("-c", "--svm-penalty",
            action="store", type="float",dest="C",default=1.),
            make_option(None, "--smlr-lm",
            action="store", type="float",dest="smlrLm",default=0.01,
            help="SMLR Lambda penalty - higher values lead to more sparse solutions"),
            make_option(None, "--smlr-convergence-tol",
            action="store", type="float",dest="smlrConvergenceTol",default=0.001,
            help="SMLR convergence tolerance - a value of 0.01 can speed up training considerably at the expense of ~2% accuracy loss"),
            make_option("-m", "--mode",
            action="store", type="choice",choices=optChoicesMode,
            dest="mode",default="train",help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            make_option("-e", "--method",
            action="store", type="choice",choices=("svm","knn","knn-svm","smlr"),dest="method",default="svm"),
            make_option("-k", "--knn-k",
            action="store", type="int",dest="knnK",default=5,help="Number of nearest neighbours in KNN algorithm"),
            make_option("-d", "--knn-max-dist",
            action="store", type="float",dest="knnMaxDist",default=None),
            make_option(None, "--kernel",
            action="store", type="choice",choices=("lin","rbf"),dest="kernel",default="lin"),
            make_option(None, "--rbf-width",
            action="store", type="float",dest="rbfWidth",default=1.),
            make_option("-l", "--train-labels",
            action="store", type="string",dest="trainLabels",default=None,
            help="Training labels to process here (all by default)"),
            make_option(None, "--balance-train-counts",
            action="store", type="int",dest="balanceTrainCounts",default=-2,
            help="Balance training samples to a given count by random subsampling. "+\
                    "if 0 - balance to the min class count, if -1 - only shuffle, if <-1 - do nothing, "+\
                    "else - to this count, if 'median' - to median count"),
            make_option(None, "--balance-random-seed",
            action="store", type="int",dest="balanceRndSeed",default=76576,
            help="Random seed to use if training set balancing is requested. This is needed because some methods (kernel svm) "+\
                    "store references to the training phase samples in the final model. This seed must be the same in the prediction/testing phase."),
            make_option(None, "--num-train-jobs",
            action="store", type="int",dest="numTrainJobs",default=1,
            help="Maximum number of training jobs to run in parallel"),
            make_option("-t", "--predict-thresh",
            action="store", type="string",dest="thresh",default="-2000"),
            make_option("-o", "--model-name",
            action="store", type="string",dest="modelRoot",default="mod"),
            make_option("-p", "--pred-file",
            action="store", type="string",dest="predFile",default=None),
            make_option("-r", "--perf-file",
            action="store", type="string",dest="perfFile",default=None),
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
        if isinstance(options.thresh,str):
            thresh = [ float(t) for t in options.thresh.split(',') ]
        else:
            thresh = options.thresh
        if len(thresh) == 3:
            thresh = numpy.arange(thresh[0],thresh[1],thresh[2])
        elif len(thresh) == 1:
            thresh = numpy.array(thresh,dtype='f4')
        else:
            parser.error("--predict-thresh value must be 'start,end,step' or 'value', received %s" % options.thresh)
        options.thresh = thresh
        if options.trainLabels is not None:
            if isinstance(options.trainLabels,str):
                options.trainLabels = n.asarray([ int(l) for l in options.trainLabels.split(',') ])
        if options.mode == "test":
            if options.perfFile is None:
                options.perfFile = "perf.pkl"
        if options.mode == "predict":
            if options.predFile is None:
                options.predFile = "pred.pkl"

    def loadLabels(self):
        opt = self.opt
        labels = opt.labels
        assert not isinstance(labels,str),"Need a sequence of label files"
        if labels is None:
            labels = [ lf+".idlab" for lf in opt.inFeat ]
        idLab = loadIdLabelsMany(fileNames=labels)
        #DEBUG:
        #idLab = idLab.balance(100)
        self.idLab = idLab
    
    def loadFeatures(self):
        opt = self.opt
        idLab = self.idLab
        assert not isinstance(opt.inFeat,str),"Need a sequence of feature files"
        assert len(opt.inFeat) > 0
        if options.debug > 0:
            timer = Timer()
            print "Loading features"
        dataFeat = loadSparseSeqsMany(opt.inFeat,idLab=idLab,format=opt.inFeatFormat)
        if options.debug > 0:
            print "Features loaded in %.3f sec" % (timer(),)
        idLabSplits = idLab.getSplits()
        feat = []
        lab = []
        data = []
        iSplit = 0
        for split in sorted(idLabSplits):
            idl = idLabSplits[split]
            if iSplit == 0: #training set only
                if opt.balanceTrainCounts >= -1:
                    assert opt.balanceRndSeed is not None, "Random seed has to be set when balancing the training data"
                    idl = idl.balance(maxCount=opt.balanceTrainCounts,rndSeed=opt.balanceRndSeed)
            splitData = idl.selData(dataFeat)
            splitFeat = convSparseToShog(data=splitData,delFeature=True)
            data.append(splitData)
            feat.append(splitFeat)
            lab.append(splitData["label"].astype('i4'))
            print "DEBUG: len(data[%s] = %s, len(feat[%s]) = %s, len(lab[%s]) = %s" % (iSplit,len(data[iSplit]),iSplit,feat[iSplit].get_num_vectors(),iSplit,len(lab[iSplit]))
            iSplit += 1
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


    def compute(self):
        feat = self.feat
        lab = self.lab
        opt = self.opt
        
        trainLabels = opt.trainLabels

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
            if opt.method == "svm":
                svmMul = SvmOneVsAll(maxLabel=maxLab)
                svmMul.setLab(lab[0])
                svmMul.setFeat(feat[0])
                svmMul.setSvmStore(modStore)
                svmMul.trainMany(trainLabels=trainLabels)
            elif opt.method == "smlr":
                import mvpa.datasets
                from mvpa.clfs.smlr import SMLR
                mv_data = mvpa.datasets.Dataset(samples=feat[0].get_full_feature_matrix().transpose(),labels=lab[0])
                clf = SMLR(lm=opt.smlrLm,convergence_tol=opt.smlrConvergenceTol)
                clf.train(mv_data)
                makedir(opt.modelRoot)
                dumpObj(clf,pjoin(opt.modelRoot,"smlr"))

        elif opt.mode in ("test","predict"):
            if opt.method == "svm":
                if trainLabels is None:
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
            elif opt.method == "smlr":
                clf = loadObj(pjoin(opt.modelRoot,"smlr"))
                labPred = n.asarray(clf.predict(feat[2].get_full_feature_matrix().transpose()),dtype='i4')
                labPred.shape = (1,len(labPred))
                return Struct(labPred=labPred,param=n.rec.fromarrays([thresh[0:1]],names="thresh"))
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

    def listAllTrainLabels(self):
        opt = self.opt
        idLab = self.idLab
        idLabSplits = idLab.getSplits()
        idLabTr = idLabSplits[sorted(idLabSplits)[0]]
        return n.asarray(sorted(idLabTr.getLabToRec()))

    def doWork(self,**kw):
        opt = self.opt
        print "Classifier options: \n%s\n" % opt
        assert not (opt.mode == "test" and opt.perfFile is None)
        assert not (opt.mode == "predict" and opt.predFile is None)
        self.loadLabels()
        if opt.trainLabels is None:
            trainLabelsAll = self.listAllTrainLabels()
            # If the method training phase has not been parallelized,
            # just train immediately in a single process
            ## @todo Parallelize test/predict modes by splitting the dataset
            if opt.method not in ["svm"]:
                if opt.mode == "trainScatter":
                    opt.numTrainJobs = 1
            if opt.mode == "trainScatter":
                # I am a dispatch job - split work into sets of labels and dispatch
                numJobs = min(opt.numTrainJobs,len(trainLabelsAll))
                if len(trainLabelsAll) < 3:
                    numJobs = 1
                trainLabelsJobs = n.array_split(trainLabelsAll,numJobs)
                # clean up all old models in this dispatch job -
                # at this point we know that all models will be built again
                rmdir(opt.modelRoot)
                jobs = []
                for trainLabels in trainLabelsJobs:
                    clOpt = copy(opt)
                    clOpt.mode = "train"
                    clOpt.trainLabels = trainLabels
                    clApp = ClassifierApp(opt=clOpt)
                    jobs.extend(clApp.run(**kw))
                return jobs
            elif opt.mode == "train":
                opt.trainLabels = trainLabelsAll
        # If we got here, we do an actual work
        self.loadFeatures()
        pred = self.compute()
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

