### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


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
        return self.crossVal(**kw)

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
        This file can be updated as computation progresses.
        It can hold file names for final result files etc."""
        return pjoin(self.cwd,"opt.pkl")

    def trainTest(self,testSplit,**kw):
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
        spClOpt.mode = "trainScatter"
        spClOpt.runMode = self.opt.runMode
        spClOpt.cwd = spDir
        spOptFile = self.getSplitOptFileName(testSplit)
        spApp = ClassifierApp(opt=spClOpt)
        spClOpt = spApp.getOpt() #get the missing options filled with defaults
        #that just saves a copy for us - App will create its own when it 
        #submits a batch job (unique and slightly modified)
        dumpObj(spClOpt,spOptFile)
        kw = kw.copy()
        jobs = spApp.run(**kw)
        spClOpt.mode = "test"
        spApp = ClassifierApp(opt=spClOpt)
        kw["depend"] = jobs
        jobs = spApp.run(**kw)
        return jobs

    def gatherSplits(self,**kw):
        """Collect results from training/testing all splits and compute total performance metrics"""
        opt = self.opt
        clOpt = self.clOpt
        labPred = []
        spPerfs = []
        param = None
        idPred = []
        if options.debug > 0:
            print "CV results:"
        for (spDir,split) in self.listSplitDirNames():
            spClOpt = loadObj(self.getSplitOptFileName(split))
            spPred = loadObj(spClOpt.predFile)
            labPred.append(spPred.labPred)
            #param (e.g. decision thresholds) is assumed to be identical for all splits
            param = spPred.param
            idPred.append(spPred.idPred)
            spPerf = loadObj(spClOpt.perfFile)
            spPerf.joinParam(n.asarray([(split,)],dtype=[("split","i4")]))
            spPerfs.append(spPerf)
        spPerfs = spPerfs[0].concatenate(spPerfs)
        labPred = n.column_stack(labPred)
        idPred = n.concatenate(idPred)
        pred = Predictions(labPred=labPred,param=param,idPred=idPred)
        opt.setdefault("predFile",pjoin(self.cwd,"pred.pkl"))
        opt.setdefault("perfFile",pjoin(self.cwd,"perf.pkl"))
        opt.setdefault("perfFileCl",pjoin(self.cwd,"perf.cl.pkl"))
        dumpObj(opt,self.getOptFileName())
        if opt.predFile is not None:
            dumpObj(pred,opt.predFile)
        perf = pred.calcPerfMetrics(idLab=self.idLab,
                confMatrFileStem=opt.perfFile if clOpt.exportConfMatr else None,
                keepConfMatr=clOpt.saveConfMatr)
        dumpObj(perf,opt.perfFile)
        dumpObj(spPerfs,opt.perfFileCl)


    def crossVal(self,**kw):
        if self.opt.mode == "scatter":
            jobs = []
            for split in sorted(self.idLab.getSplitToRec()):
                jobs.extend(self.trainTest(testSplit=split,**kw))
            gtOpt = copy(self.opt)
            gtOpt.mode = "gather"
            gtApp = self.factory(opt=gtOpt)
            kw = kw.copy()
            kw["cwd"] = self.cwd
            kw["depend"] = jobs
            jobs = gtApp.run(**kw)
            return jobs
        elif self.opt.mode == "gather":
            self.gatherSplits(**kw)



#if __name__ == "__main__":
#    runAppAsScript(CrossValidatorApp)

