from MGT.App import *
from MGT.ClassifierApp import ClassifierApp
from MGT.Svm import *
from MGT.PredProcessor import Predictions

class CrossValidatorApp(App):
    """Class that performs cross validation"""

    batchDepModes = ("scatter",)

    def doWork(self,**kw):
        opt = self.opt
        self.clOpt = opt.clOpt
        self.cwd = opt.get("cwd",os.getcwd())
        self.loadIdLabels()
        return self.crossVal()

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
        return pjoin(self.getSplitDirName(split),"opt.pkl")

    def getOptFileName(self):
        """Return stable file name for options dump.
        This file can be updated as compputation progresses.
        It can hold file names for final result files etc."""
        return pjoin(self.cwd,"opt.pkl")

    def trainTest(self,testSplit):
        idLab = self.idLab
        spDir = self.getSplitDirName(testSplit)
        makedir(spDir)
        idLabRec = idLab.getRecords()
        spIdLab = idLabRec.copy()
        spIdLab["split"] = 1
        spIdLab["split"][idLabRec["split"]==testSplit] = 2
        spIdLabFile = pjoin(spDir,"idlab.pkl")
        saveIdLabelRecords(records=spIdLab,fileName=spIdLabFile)
        spClOpt = copy(self.clOpt)
        spClOpt.labels = [spIdLabFile]
        spClOpt.predFile = pjoin(spDir,"pred.pkl")
        spClOpt.perfFile = pjoin(spDir,"perf.pkl")
        spClOpt.mode = "train"
        spClOpt.runMode = self.opt.runMode
        spOptFile = self.getSplitOptFileName(testSplit)
        #that just saves a copy for us - App will create its own when it 
        #submits a batch job (unique and slightly modified)
        dumpObj(spClOpt,spOptFile)
        spApp = ClassifierApp(opt=spClOpt)
        jobs = spApp.run(cwd=spDir)
        spClOpt.mode = "test"
        spApp = ClassifierApp(opt=spClOpt)
        jobs = spApp.run(cwd=spDir,depend=jobs)
        return jobs

    def gatherSplits(self):
        """Collect results from training/testing all splits and compute total performance metrics"""
        opt = self.opt
        clOpt = self.clOpt
        labPred = []
        param = None
        idPred = []
        for (spDir,split) in self.listSplitDirNames():
            spClOpt = loadObj(self.getSplitOptFileName(split))
            spPred = loadObj(spClOpt.predFile)
            labPred.append(spPred.labPred)
            #param (e.g. decision thresholds) is assumed to be identical for all splits
            param = spPred.param
            idPred.append(spPred.idPred)
        labPred = n.column_stack(labPred)
        idPred = n.concatenate(idPred)
        pred = Predictions(labPred=labPred,param=param,idPred=idPred)
        opt.setdefault("predFile",pjoin(self.cwd,"pred.pkl"))
        opt.setdefault("perfFile",pjoin(self.cwd,"perf.pkl"))
        dumpObj(opt,self.getOptFileName())
        if opt.predFile is not None:
            dumpObj(pred,opt.predFile)
        perf = pred.calcPerfMetrics(idLab=self.idLab,
                confMatrFileStem=opt.perfFile if clOpt.exportConfMatr else None,
                keepConfMatr=clOpt.saveConfMatr)
        dumpObj(perf,opt.perfFile)


    def crossVal(self):
        if self.opt.mode == "scatter":
            jobs = []
            for split in sorted(self.idLab.getSplitToRec()):
                jobs.extend(self.trainTest(testSplit=split))
            gtOpt = copy(self.opt)
            gtOpt.mode = "gather"
            gtApp = self.factory(opt=gtOpt)
            jobs = gtApp.run(cwd=self.cwd,depend=jobs)
            return jobs
        elif self.opt.mode == "gather":
            self.gatherSplits()



#if __name__ == "__main__":
#    runAppAsScript(CrossValidatorApp)

