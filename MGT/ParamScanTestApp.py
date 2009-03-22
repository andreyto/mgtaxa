from MGT.Taxa import *
from MGT.SeqIO import pullNCBISeqIds
from MGT.Svm import *
from MGT.App import *
from MGT.ClassifierApp import *
from MGT.ParamScanApp import *
from MGT.SeqFeaturesApp import *

from MGT.DirStore import *

class ParamScanTestApp(App):
    """App-derived class for one-stop feature creation - parameter scan - final testing"""

    batchDepModes = ("parscan","test")

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.parScanName = "2"
        opt.featName = "kmer_8" #"kmer_8a" #"kmerlad_4f_nnorm_rndGenGcFixed" #"kmerlad_4f" #"wd_2_10"
        clOpt = opt.setdefault("clOpt",Struct())
        psOpt = opt.setdefault("psOpt",Struct())
        ftOpt = opt.setdefault("ftOpt",Struct())
        clOpt.thresh = n.arange(-2.,1.01,0.1)
        clOpt.C = 4.
        #clOpt.saveConfMatr = True
        #clOpt.exportConfMatr = True
        clOpt.method = "svm" #"svm" "knn"
        clOpt.kernel = "lin" #"rbf"
        clOpt.knnK = 10
        ClassifierApp.fillWithDefaultOptions(clOpt)
        pgen = ParamGridGen()
        #psOpt.params = pgen.add("C",[1]).grid() #(-10,10,2)(0,1,2)
        #psOpt.params = pgen.add("C",pgen.p2(6,16,2)).add("rbfWidth",pgen.p2(-4,10,2)).grid()
        psOpt.params = pgen.add("C",pgen.p2(-8,10,2)).grid() #(-10,10,2)(0,1,2)
        #psOpt.params = pgen.add("knnK",pgen.lin(1,20,2)).grid()
        #psOpt.params = pgen.add("knnMaxDist",pgen.lin(0.01,0.1,0.01)).grid()
        psOpt.clOpt = clOpt
        #ftOpt.inSeq = "samp.rndGenGcFixed"
        ftOpt.featType = "kmer" #"wdh" "kmerlad" #"kmer"
        if ftOpt.featType == "wdh":
            ftOpt.revCompl = "merge" #"addcol" #"forward"
        elif ftOpt.featType == "kmerlad":
            ftOpt.kmerLen = 9
            ftOpt.revCompl = "addcol" #"forward"
            ftOpt.norm = NORM_POLICY.EXPECT | NORM_POLICY.EU_ROW
        elif ftOpt.featType == "kmer":
            ftOpt.kmerLen = 8
            ftOpt.revCompl = "merge" #"forward"
            ftOpt.norm = NORM_POLICY.FREQ | NORM_POLICY.EU_ROW
        ftOpt.balance = -2 #do not even shuffle (somehow it takes a while) - we do it anyway when making idlabs
    
    def init(self):
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
    
    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "feat":
            self.prepFeat()
        elif opt.mode == "featPred":
            self.prepFeatPred()
        elif opt.mode == "idlabs":
            self.makeIdLabs()
        elif opt.mode == "parscan":
            return self.parScan(parScanName=opt.get("parScanName",1))
        elif opt.mode == "test":
            return self.test()
        elif opt.mode == "predict":
            self.predict()
        else:
            raise ValueError(opt.mode)


    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def getFeatStore(self,featName):
        return self.store.loadStore(name=featName)

    def prepFeat(self):
        ftOpt = self.opt.ftOpt
        featStore = self.store.featStore(name=self.opt.featName,opt=ftOpt,mode="w")
        featStore.fromOther(sampStore=self.store)

    def makeIdLabs(self,idLabsName=None):
        opt = self.opt
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)
        idLabsSamp = store.loadIdLabs()
        splitSet = idLabsSamp.splitSet()
        idLabsSampCv = idLabsSamp.selBySplits(splitSet-set([opt.splitIdTest,opt.splitIdTrainPred]))
        idLabsSampCv = idLabsSampCv.selAcrossSplits()
        idLabsSampCv = idLabsSampCv.balance(maxCount="median")
        store.saveIdLabs(idLabsSampCv,name=nameIL+'.cv')
        print "idLabsSampCv = ",idLabsSampCv.count()
        featStore = self.getFeatStore(opt.featName)
        idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsSampCv,name=nameIL+".cv",action="idfilt")
        #idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsFeat,name=nameIL+".cv",action="balance")
        idLabsRec = idLabsSamp.getRecords().copy()
        splitRec = idLabsRec["split"]
        splitRec[splitRec == opt.splitIdTest] = -1
        splitRec[splitRec != -1] = 1
        splitRec[splitRec == -1] = 2
        idLabsSampTest = IdLabels(records=idLabsRec)
        ## @todo The following is a cludge until we figure out a better perfomance metric,
        ## otherwise training labels w/o any testing samples will get zero specificities even
        ## after a single false positive and artificially drag down the average specificity.
        idLabsSampTest = idLabsSampTest.selAcrossSplits()
        idLabsSampTest = idLabsSampTest.balance(maxCount="median")
        store.saveIdLabs(idLabsSampTest,name=nameIL+'.ts')
        print "idLabsSampTest = ",idLabsSampTest.count()
        idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsSampTest,name=nameIL+".ts",action="idfilt")
        #idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsFeat,name=nameIL,action="balance")

    def parScan(self,parScanName=None,idLabsName=None):
        opt = self.opt
        psOpt = opt.psOpt
        psOpt.runMode = opt.runMode
        jobs = []
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)+'.cv'
        featStore = self.getFeatStore(opt.featName)
        jobs.extend(featStore.parScan(opt=psOpt,name=parScanName,idLabsName=nameIL))
        return jobs


    def test(self,testName=None,idLabsName=None):
        opt = self.opt
        clOpt = opt.psOpt.clOpt
        clOpt.runMode = opt.runMode
        jobs = []
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)+'.ts'
        featStore = self.getFeatStore(opt.featName)
        jobs.extend(featStore.test(opt=clOpt,name=testName,idLabsName=nameIL))
        return jobs


    def prepFeatPred(self):
        opt = self.opt
        ftOpt = self.opt.ftOpt
        sampStorePred = SampStore.open(path=opt.sampStorePred)
        featStore = sampStorePred.featStore(name=self.opt.featName,opt=ftOpt,mode="w")
        featStore.fromOther(sampStore=sampStorePred)
        idLabs = IdLabels.fromIdFile(featFile=featStore.getSampFilePath(),label=0,split=3)
        featStore.saveIdLabs(idLabs)

    def predict(self):
        from MGT.Taxa import loadTaxaTree
        opt = self.opt
        store = self.store
        clOpt = opt.psOpt.clOpt
        clOpt.thresh = [0.] #[-2000]
        clOpt.runMode = opt.runMode
        featStoreTrain = self.getFeatStore(opt.featName)
        clOpt.modelRoot = pjoin(featStoreTrain.getPath(),"test/1",clOpt.modelRoot)
        sampStorePred = SampStore.open(path=opt.sampStorePred)
        featStore = sampStorePred.loadStore(name=opt.featName)
        jobs = []
        jobs.extend(featStore.predict(opt=clOpt))
        idLab = store.loadIdLabs()
        labToName = idLab.getLabToName()
        pred = featStore.predictStore().loadObj("pred")
        idToName = {}
        for id,lab in zip(pred.idPred["id"],pred.labPred[0]):
            idToName[id] = labToName[lab]
        taxaTree = loadTaxaTree()
        idToLin = {}
        for id,labn in idToName.items():
            if not isinstance(labn,str):
                idToLin[id] = taxaTree.getNode(labn).lineageStr()
            else:
                idToLin[id] = labn
        out = open(featStore.predictStore().getFilePath("idlin.csv"),'w')
        for (id,lin) in sorted(idToLin.items()):
            out.write("%s\t%s\n" % (id,lin))
        out.close()
        linCnt = "\n".join(["%s\t%s" % y for y in sorted(binCount(idToLin.values()).items(),key=lambda x:-x[1])])
        out = open(featStore.predictStore().getFilePath("lincnt.csv"),'w')
        out.write(linCnt)
        out.write("\n")
        out.close()
        return jobs
