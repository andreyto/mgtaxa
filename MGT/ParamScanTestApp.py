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

    batchDepModes = ("parscan","test","train")

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        pass
    
    def init(self):
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
    
    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "feat":
            return self.prepFeat(**kw)
        elif opt.mode == "featPred":
            return self.prepFeatPred(**kw)
        elif opt.mode == "idlabs":
            return self.makeIdLabs(**kw)
        elif opt.mode == "parscan":
            return self.parScan(parScanName=opt.get("parScanName",1),**kw)
        elif opt.mode == "test":
            return self.test(**kw)
        elif opt.mode == "train":
            return self.train(**kw)
        elif opt.mode == "predict":
            return self.predict(**kw)
        else:
            raise ValueError(opt.mode)


    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def getFeatStore(self,featName):
        return self.store.loadStore(name=featName)

    def prepFeat(self,**kw):
        ftOpt = self.opt.ftOpt
        featStore = self.store.featStore(name=self.opt.featName,opt=ftOpt,mode="w")
        return featStore.fromOther(sampStore=self.store,**kw)

    def makeIdLabs(self,idLabsName=None,**kw):
        opt = self.opt
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)
        idLabsSamp = store.loadIdLabs()
        splitSet = idLabsSamp.splitSet()
        idLabsSampCv = idLabsSamp.selBySplits(splitSet-set([opt.splitIdTest,opt.splitIdTrainPred]))
        idLabsSampCv = idLabsSampCv.renumSplits()
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
        idLabsSampTestDbg = idLabsSampTest.selByLabs(idLabsSampTest.labSet()-idLabsSampCv.labSet())
        #idLabsSampTestDbg = idLabsSampTest.selByLabs(idLabsSampCv.labSet())
        store.saveIdLabs(idLabsSampTestDbg,name=nameIL+'.dbg'+'.ts')
        print "idLabsSampTestDbg = ",idLabsSampTestDbg.count()
        idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsSampTestDbg,name=nameIL+'.dbg'+".ts",action="idfilt")
        # Make prediction training IdLabels - everything goes to a single split 1 and then balanced
        idLabsRec = idLabsSamp.getRecords().copy()
        idLabsRec["split"] = 1
        idLabsSampTrain = IdLabels(records=idLabsRec)
        idLabsSampTrain = idLabsSampTrain.balance(maxCount="median")
        store.saveIdLabs(idLabsSampTrain,name=nameIL+'.tr')
        print "idLabsSampTrain = ",idLabsSampTrain.count()
        idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsSampTrain,name=nameIL+".tr",action="idfilt")
        #idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsFeat,name=nameIL,action="balance")

    def parScan(self,parScanName=None,idLabsName=None,**kw):
        opt = self.opt
        psOpt = opt.psOpt
        psOpt.runMode = opt.runMode
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)+'.cv'
        featStore = self.getFeatStore(opt.featName)
        return featStore.parScan(opt=psOpt,name=parScanName,idLabsName=nameIL,**kw)


    def test(self,testName=None,idLabsName=None,**kw):
        opt = self.opt
        clOpt = opt.tsOpt.clOpt.copy()
        clOpt.runMode = opt.runMode
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)+'.dbg.ts' #TMP:'.dbg.ts'
        featStore = self.getFeatStore(opt.featName)
        return featStore.test(opt=clOpt,name=testName,idLabsName=nameIL,**kw)

    def train(self,testName=None,idLabsName=None,**kw):
        """Train on all available samples to use models for prediction of future unknown samples"""
        opt = self.opt
        clOpt = opt.prOpt.clOpt
        clOpt.runMode = opt.runMode
        store = self.store
        nameIL = store.getIdLabsName(idLabsName)+'.tr'
        featStore = self.getFeatStore(opt.featName)
        return featStore.train(opt=clOpt,name=testName,idLabsName=nameIL,**kw)

    def prepFeatPred(self,**kw):
        opt = self.opt
        ftOpt = self.opt.ftOpt.copy()
        ftOpt.runMode = "inproc"
        sampStorePred = SampStore.open(path=opt.sampStorePred)
        featStore = sampStorePred.featStore(name=self.opt.featName,opt=ftOpt,mode="w")
        featStore.fromOther(sampStore=sampStorePred)
        idLabs = IdLabels.fromIdFile(featFile=featStore.getSampFilePath(),label=0,split=3)
        featStore.saveIdLabs(idLabs)

    def predict(self,**kw):
        from MGT.Taxa import loadTaxaTree
        opt = self.opt
        store = self.store
        clOpt = opt.prOpt.clOpt.copy()
        assert len(clOpt.thresh) == 1
        clOpt.runMode = "inproc" #we would need to split this into two methods for "batchDep"
        featStoreTrain = self.getFeatStore(opt.featName)
        clOpt.modelRoot = pjoin(featStoreTrain.getPath(),"train/1",clOpt.modelRoot)
        sampStorePred = SampStore.open(path=opt.sampStorePred)
        featStore = sampStorePred.loadStore(name=opt.featName)
        featStore.predict(opt=clOpt,name=opt.prOpt.name)
        idLab = store.loadIdLabs()
        labToName = idLab.getLabToName()
        predStore = featStore.predictStore(name=opt.prOpt.name)
        pred = predStore.loadObj("pred")
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
        out = open(predStore.getFilePath("idlin.csv"),'w')
        for (id,lin) in sorted(idToLin.items()):
            out.write("%s\t%s\n" % (id,lin))
        out.close()
        linCnt = "\n".join(["%s\t%s" % y for y in sorted(binCount(idToLin.values()).items(),key=lambda x:-x[1])])
        out = open(predStore.getFilePath("lincnt.csv"),'w')
        out.write(linCnt)
        out.write("\n")
        out.close()

