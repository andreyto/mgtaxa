### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Read a file with sequence features and create a file with sparse real features based on word distance histogram."""
from MGT.Taxa import *
from MGT.Svm import *
from MGT.FastaIO import *
from MGT.Shogun.Util import *
from shogun.Features import *
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-t", "--tax-ids",
        action="store", type="string",dest="taxids"),
        make_option("-o", "--out-file",
        action="store", type="string",dest="outFile"),
        make_option("-s", "--with-subnodes",
        action="store_true", dest="withSubnodes",default=True),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

taxids = [ int(x.strip()) for x in opt.taxids.split(',') ]

mgtDbDir = "/home/atovtchi/work/mgtdata"

taxaDir=os.path.join(mgtDbDir,"taxonomy.new")

giToTaxa = loadGiTaxBin(os.path.join(taxaDir,"gi_taxid.pkl.gz"))

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))

downloadTaxaGenbankSeqs(taxids=taxids,outFile=opt.outFile,
        taxaTree=taxaTree,giToTaxa=giToTaxa,withSubnodes=opt.withSubnodes)

