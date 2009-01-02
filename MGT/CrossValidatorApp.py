from MGT.App import *

class CrossValidator(App):
    """Class that performs cross validation"""

    def run(self,**kw):
        opt = self.opt
        if opt.batch:
            return self.runBatch(**kw)
        self.clOpt = opt.clOpt
        self.cwd = opt.get("cwd",os.getcwd())
        self.loadIdLabels()
        self.crossVal()

    def loadIdLabels(self):
        self.idLab = loadIdLabelsMany(fileNames=self.clOpt.labels)

    def getSplitDirName(self,split):
        return pjoin(self.cwd,"split-%00i" % split)

    def listSplitDirNames(self):
        from glob import glob
        import re
        rex = re.compile(r"split-([0-9]+)/*$")
        return [ (d,int(rex.search(d).group(1))) for d in glob(pjoin(self.cwd,"split-*")) ]

    def getSplitOptFileName(self,split):
        return pjoin(self.getSplitDirName(),"opt.pkl")

    def getOptFileName(self):
        """Return stable file name for options dump.
        This file can be updated as compputation progresses.
        It can hold file names for final result files etc."""
        return pjoin(self.cwd,"opt.pkl")

    def trainTest(self,testSplit):
        idLab = self.idLab
        spDir = self.getSplitDirName(testSplit)
        makedir(spDir)
        spIdLab = idLab.getRecords().copy()
        spIdLab["split"] = 1
        spIdLab["split"][idLab["split"]==testSplit] = 2
        spIdLabFile = pjoin(spDir,"idlab.pkl")
        saveIdLabelRecords(records=spIdLab,fileName=psIdLabFile)
        spClOpt = copy(self.clOpt)
        spClOpt.labels = [spIdLabFile]
        spClOpt.predFile = pjoin(spDir,"pred.pkl")
        spClOpt.perfFile = pjoin(spDir,"perf.pkl")
        spClOpt.mode = "train"
        spOptFile = self.getSplitOptFileName(testSplit)
        #that just saves a copy for us - App will create its own when it 
        #submits a batch job (unique and slightly modified)
        dumpObj(spClOpt,spOptFile)
        spApp = ClassifierApp(**spClOpt)
        jobs = spApp.runBatch(cwd=spDir)
        spClOpt.mode = "test"
        spApp = ClassifierApp(**spClOpt)
        jobs = spApp.runBatch(cwd=spDir,depend=jobs)
        return jobs

    def gatherTest(self):
        """Collect results from training/testing all splits"""
        opt = self.opt
        pred = Struct(pred=Struct(labPred=[],param=None),idPred=[]) 
        for (spDir,split) in self.listSplitDirNames():
            spClOpt = loadObj(self.getSplitOptFileName(split))
            spPred = loadObj(spClOpt.predFile)
            pred.pred.labPred.append(spPred.pred.labPred)
            #param (e.g. decision thresholds) is assumed to be identical for all splits
            pred.pred.param = spPred.pred.param
            pred.idPred.append(spPred.idPred)
        pred.pred.labPred = n.column_stack(pred.pred.labPred)
        pred.idPred = n.column_stack(pred.idPred)
        opt.setdefault("predFile",pjoin(self.cwd,"pred.pkl"))
        opt.setdefault("perfFile",pjoin(self.cwd,"perf.pkl"))
        dumpObj(pred,opt.predFile)
        dumpObj(opt,self.getOptFileName())





    def crossVal(self):
        jobs = []
        for split in sorted(self.idLab.getSplitToRec()):
            jobs.extend(self.trainTest(testSplit=split))
        gtOpt = copy(self.opt)
        gtOpt.mode = "gatherTest"
        gtApp = self.getFactory()(**gtOpt)
        jobs = gtApp.runBatch(cwd=self.cwd,depend=jobs)
        return jobs



if __name__ == "__main__":
    app = CrossValidatorApp()
    app.run()

