"""Split a file with sparse real features"""

from MGT.Shogun import *
from MGT.Shogun.Util import *
from shogun.Features import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-feat",
        action="store", type="string",dest="inFeat"),
        make_option("-o", "--out-feat",
        action="store", type="string",dest="outFeat"),
        make_option("-r", "--ratio",
        action="store", type="float",dest="splitRatio",default=0.5),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

feat = SparseRealFeatures()
lab = feat.load_svmlight_file(opt.inFeat)
print "Loaded %d samples" % feat.get_num_vectors()
splitTwo(feat,lab,opt.splitRatio,opt.outFeat)

