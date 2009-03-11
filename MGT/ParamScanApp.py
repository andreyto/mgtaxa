from MGT.App import *
from MGT.ClassifierApp import ClassifierApp
from MGT.CrossValidatorApp import CrossValidatorApp
from MGT.Svm import *
from MGT.PredProcessor import Predictions

__all__ = ["ParamGridGen","ParamScanApp"] 

def scalePow(base,start,stop,step=1):
    """Return an array that represents sequential powers of a given base"""
    return n.power(float(base),n.arange(start,stop,step))

def scalePow2(start,stop,step=1):
    """Return an array that represents sequential powers of two"""
    return scalePow(2.,start,stop,step)

def scaleLin(start,stop,step=1):
    """Return a linear scale array"""
    return n.arange(start,stop,step)

class ParamGridGen:
    """Class that allows to add N 1D scales and then obtain a 1D record array enumerating all N-D values.
    Example of use:
    g = ParamGridGen()
    g.add("C",g.p2(-6,4,2)).add("gamma",g.lin(0.01,10,0.2)).grid()"""

    ## Short aliases to scale generating functions
    p2 = staticmethod(scalePow2)
    lin = staticmethod(scaleLin)

    def __init__(self):
        self.scales = []
    
    def add(self,name,values):
        """Add a scale (dimention) with a given name and values.
        Use scaleXXX set of functions to construct the values."""
        self.scales.append(Struct(name=name,values=n.asarray(values)))
        return self

    def grid(self):
        """Return a 1D record array enumerating all combinations of scales' values"""
        scales = self.scales
        dt = [ (s.name,s.values.dtype) for s in scales ]
        vals = [ s.values for s in scales ]
        sizes = [ len(val) for val in vals ]
        nel = n.product(sizes)
        out = n.zeros(nel,dtype=dt)
        for (ind,i) in it.izip(n.ndindex(*sizes),it.count()):
            out[i] = tuple([ v[j] for (v,j) in it.izip(vals,ind) ])
        return out



class ParamScanApp(App):
    """Class that runs cross-validation for different combinations of classifier's hyper-parameters."""

    batchDepModes = ("scatter",)

    def doWork(self,**kw):
        opt = self.opt
        self.clOpt = opt.clOpt
        self.cwd = opt.get("cwd",os.getcwd())
        if opt.mode == "scatter":
            jobs = []
            for (param,idParam) in it.izip(opt.params,it.count()):
                jobs.extend(self.crossVal(param=param,idParam=idParam))
            gtOpt = copy(opt)
            gtOpt.mode = "gather"
            gtApp = self.factory(opt=gtOpt)
            jobs = gtApp.run(cwd=self.cwd,depend=jobs)
            return jobs
        elif opt.mode == "gather":
            self.gather()

    def getCvDirName(self,idParam):
        return pjoin(self.cwd,"cv-%000i" % idParam)

    def getOptFileName(self):
        """Return stable file name for options dump.
        This file can be updated as computation progresses.
        It can hold file names for final result files etc."""
        return pjoin(self.cwd,"opt.pkl")

    def getCvOptFileName(self,idParam):
        return pjoin(self.getCvDirName(idParam),"opt.pkl")
    
    def crossVal(self,param,idParam):
        opt = self.opt
        cvDir = self.getCvDirName(idParam)
        makedir(cvDir)
        cvClOpt = copy(self.clOpt)
        for parName in param.dtype.names:
            parVal = param[parName]
            setattr(cvClOpt,parName,parVal.item()) #convert numpy scalar to python built-in type
        cvOpt = Struct(clOpt=cvClOpt)
        cvOpt.predFile = pjoin(cvDir,"pred.pkl")
        cvOpt.perfFile = pjoin(cvDir,"perf.pkl")
        cvOpt.mode = "scatter"
        cvOpt.runMode = opt.runMode
        cvApp = CrossValidatorApp(opt=cvOpt)
        jobs = cvApp.run(cwd=cvDir)
        return jobs

    def gather(self):
        """Collect results from cross-validation jobs and save performance metrics in a single PerfMetricsSet"""
        opt = self.opt
        clOpt = self.clOpt
        perfAll = []
        perfClAll = []
        for (param,idParam) in it.izip(opt.params,it.count()):
            cvDir = self.getCvDirName(idParam)
            cvOpt = loadObj(self.getCvOptFileName(idParam))
            perf = loadObj(cvOpt.perfFile)
            perf.joinParam(param)
            perfAll.append(perf)
            perfCl = loadObj(cvOpt.perfFileCl)
            perfCl.joinParam(param)
            perfClAll.append(perfCl)
            
        perfAll = perfAll[0].concatenate(perfAll)
        opt.setdefault("perfFile",pjoin(self.cwd,"perf.pkl"))
        dumpObj(perfAll,opt.perfFile)
        perfAll.exportMetricsCsv(names=("senMin","speMin","senMean","speMean","acc"),out=pjoin(self.cwd,"perf.csv"))
        perfClAll = perfClAll[0].concatenate(perfClAll)
        opt.setdefault("perfFileCl",pjoin(self.cwd,"perf.cl.pkl"))
        dumpObj(perfClAll,opt.perfFileCl)
        perfClAll.exportMetricsCsv(names=("senMin","speMin","senMean","speMean","acc"),out=pjoin(self.cwd,"perf.cl.csv"))
        dumpObj(opt,self.getOptFileName())



#if __name__ == "__main__":
#    runAppAsScript(ParamScanApp)

