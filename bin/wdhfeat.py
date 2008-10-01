"""Read a file with sequence features and create a file with sparse real features based on word distance histogram."""

from MGT.Shogun.Util import *
from shogun.Features import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="store", type="string",dest="inSeq"),
        make_option("-o", "--out-feat",
        action="store", type="string",dest="outFeat"),
        make_option("-s", "--sigma",
        action="store", type="float",dest="sigma",default=200),
        make_option("-k", "--kmer-len",
        action="store", type="int",dest="kmerLen",default=2),
        make_option("-d", "--kmer-max-dist",
        action="store", type="int",dest="maxDist",default=10),
        make_option("-l", "--kmer-min-dist",
        action="store", type="int",dest="minDist",default=-1),
        make_option("-b", "--balance",
        action="store", type="int",dest="balance",default=-1),
        make_option("-u", "--other-group",
        action="store", type="int",dest="otherGroupLab",default=0),
        make_option("-r", "--rev-compl",
        action="store", type="choice",choices=("merge","forward","addcol","addrow","reverse"),
        dest="revCompl",default="merge"),
        make_option("-a", "--alphabet",
        action="store", type="choice",choices=("dna","protein"),
        dest="alphabet",default="dna"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


def otherVsRest(data,otherLabel):
    data['label'][:] = numpy.select([data['label'] == otherLabel],[1],default=2)
    return data

opt,args = getProgOptions()

data = loadSeqs(opt.inSeq)

##data = otherVsRest(data,2)

if opt.balance >= 0:
    data = balance(data,opt.balance,labTargets={opt.otherGroupLab:-1})
##TMP:
##data = splitStringFeat(data,750)
if opt.alphabet == 'dna':
    shogAlpha = DNA
    applyToFeatData(data,transDegen)
else:
    shogAlpha = PROTEIN
    assert opt.revCompl == 'forward'

rcPolicy = WH_RC_FORWARD
if opt.revCompl == "merge":
    rcPolicy=WH_RC_MERGE
elif opt.revCompl == "addcol":
    data = addRevComplCols(data)
elif opt.revCompl == "forward":
    pass
elif opt.revCompl == "addrow":
    data = addRevComplRows(data)
elif opt.revCompl == "reverse":
    data = applyRevCompl(data)
else:
    raise ValueError("Value %s for revCompl is not supported" % opt.revCompl)

print "Program options are:\n%s\n" % (opt,)

print "Computing features: len(data) = %s counts(data) = %s" % (len(data), numpy.bincount(data['label'].astype(int)))

feat_char=StringCharFeatures(shogAlpha)
feat_char.set_string_features(data['feature'].tolist())

feat_whd=WordHistogramFeatures()
feat_whd.obtain_from_char(feat_char,opt.kmerLen,opt.maxDist,opt.sigma,opt.minDist,rcPolicy)
print "WHD number of feature elements: %d" % feat_whd.get_num_elements()
feat_sparse = feat_whd.get_sparse_real_features()
print "WHD number feature dimensionality: %d" % feat_whd.get_num_features()

lab = Labels(data['label'])
feat_sparse.write_svmlight_file(opt.outFeat,lab)

svmSaveId(data['id'],opt.outFeat+'.id')

