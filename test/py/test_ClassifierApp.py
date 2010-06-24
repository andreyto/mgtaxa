### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.ClassifierApp import ClassifierApp

def makeClOpt(testDataDir):
    o,a = ClassifierApp.defaultOptions()
    o.method = "svm" #"svm" #"smlr"
    o.C = 0.001
    o.kernel = "rbf" #"rbf" "lin"
    o.inFeatFormat = "pkl" #"txt" #"pkl"
    if o.inFeatFormat == "txt":
        ext = ".libsvm"
    elif o.inFeatFormat == "pkl":
        ext = ".pkl"
    else:
        raise ValueError(o.inFeatFormat)
    origFeat = [pjoin(testDataDir,"usps-train"+ext),pjoin(testDataDir,"usps-test"+ext)]
    o.inFeat = origFeat
    #o.inFeat = [pjoin(testDataDir,"usps.libsvm.pca")]
    o.labels = [ lf+".idlab" for lf in origFeat ]
    o.balanceTrainCounts=-1
    return o

testDataDir = pjoin(options.testDataDir,"usps")

workDir = pjoin(os.getcwd(),"test_ClassifierApp.tmpdir")
rmdir(workDir)
makedir(workDir)
os.chdir(workDir)

clOpt = makeClOpt(testDataDir)
clOpt.numTrainJobs = 1 #6
clOpt.MEM = 500
clOpt.LENGTH = "fast"
clOpt.mode = "trainScatter"
clOpt.runMode = "inproc" #"inproc" #"batchDep"
print clOpt
app = ClassifierApp(opt=clOpt)
jobs = app.run(dryRun=False)
clOpt.mode = "test"
clOpt.perfFile = "perf.pkl"
app = ClassifierApp(opt=clOpt)
jobs = app.run(dryRun=False,depend=jobs)

