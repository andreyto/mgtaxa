### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.CrossValidatorApp import CrossValidatorApp
from MGT.ClassifierApp import ClassifierApp

def makeClOpt():
    o,a = ClassifierApp.defaultOptions()
    o.C = 0.03
    o.inFeat = [os.path.abspath("usps-train.libsvm"),os.path.abspath("usps-test.libsvm")]
    o.labels = [ lf+".idlab" for lf in o.inFeat ]
    return o

clOpt = makeClOpt()
print clOpt
o = Options()
o.clOpt = clOpt
o.mode = "scatter"
o.runMode = "batch"

app = CrossValidatorApp(**o.asDict())
app.run()

