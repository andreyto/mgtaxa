### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.GHSOM import GHSOM
from MGT.SOMCode import SOMCode
#from MGT.Shogun.Util import *
from shogun.Features import *
from MGT.Svm import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-feat",
        action="store", type="string",dest="inFeat"),
        make_option("-n", "--name",
        action="store", type="string",dest="name"),
        make_option("-f", "--out-format",
        action="store", type="choice",choices=("ghsom","csv","somcode"),dest="outFormat"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

feat = SparseRealFeatures()
lab = feat.load_svmlight_file(opt.inFeat).get_labels().astype('i4')
id = svmLoadId(opt.inFeat+'.id')

if opt.outFormat == "ghsom":
    som = GHSOM(name=opt.name)
    som.writeInputFromShogunSparse(feat=feat,lab=id)
if opt.outFormat == "somcode":
    som = SOMCode(name=opt.name)
    som.writeInputFromShogunSparse(feat=feat,lab=id)
elif opt.outFormat == "csv":
    mat = feat.get_full_feature_matrix().transpose()
    writer = SvmDenseFeatureWriterCsv(opt.name,writeHeader=True,nFeat=mat.shape[1])
    for iRow in xrange(mat.shape[0]):
        writer.write(lab[iRow],mat[iRow],id[iRow])
    writer.close()

