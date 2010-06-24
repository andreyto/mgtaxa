### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.ParamScanApp import ParamScanApp,ParamGridGen
from MGT.ClassifierApp import ClassifierApp

def makeClOpt(testDataDir):
    o,a = ClassifierApp.defaultOptions()
    o.C = 0.03
    o.kernel = "lin" #"rbf"
    origFeat = [pjoin(testDataDir,"usps-train.pkl"),pjoin(testDataDir,"usps-test.pkl")]
    o.inFeat = origFeat
    #o.inFeat = [pjoin(testDataDir,"usps.libsvm.pca")]
    o.labels = [ lf+".idlab" for lf in origFeat ]
    return o

testDataDir = pjoin(options.testDataDir,"usps")

workDir = pjoin(os.getcwd(),"test_ParamScanApp.tmpdir")
rmdir(workDir)
makedir(workDir)
os.chdir(workDir)

clOpt = makeClOpt(testDataDir)
clOpt.thresh = n.arange(-2.,1.,0.1)
o = Options()
o.clOpt = clOpt
o.mode = "scatter" #"gather" #"scatter"
o.runMode = "batchDep" #"batchDep" #"inproc" #"batchDep"
pgen = ParamGridGen()
#params = pgen.add("C",pgen.p2(-2,6,2)).add("rbfWidth",pgen.p2(3,16,3)).grid()
#params = pgen.add("C",pgen.p2(-2,6,2)).grid()
params = pgen.add("C",pgen.p2(-2,0,2)).grid()
o.params = params
print o
app = ParamScanApp(opt=o)
app.run()

