"""Re-balance a feature file."""

from MGT.Shogun.Util import *
from shogun.Features import *
from MGT.Svm import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-feat",
        action="store", type="string",dest="inFeat"),
        make_option("-o", "--out-feat",
        action="store", type="string",dest="outFeat"),
        make_option("-b", "--balance",
        action="store", type="int",dest="balance",default=-1),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

data = loadSeqs(opt.inFeat)

data = balance(data,opt.balance)

print "Program options are:\n%s\n" % (opt,)

out = SvmStringFeatureWriterTxt(opt.outFeat)

for rec in data:
    out.write(rec['label'],rec['feature'],rec['id'])

out.close()

