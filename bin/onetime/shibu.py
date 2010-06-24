### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Scriplets for taxonomy search for Shibu."""
from MGT.Taxa import *
from MGT.Common import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-file",
        action="store", type="string",dest="inFile"),
        make_option("-o", "--out-file",
        action="store", type="string",dest="outFile"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

mgtDbDir = "/home/atovtchi/work/mgtdata"

taxaDir=os.path.join(mgtDbDir,"taxonomy.new")

print "Loading taxonomy DB..."

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))

giToTaxa = loadGiTaxBin(os.path.join(taxaDir,"gi_taxid.pkl.gz"))

topNodes = [ taxaTree.getNode(taxid) for taxid in skingdomTaxids ]

inp = openCompressed(opt.inFile,'r')
out = openCompressed(opt.outFile,'w')

mode = "lineage"

print "Processing input data..."

iRec=0
for line in inp:
    line = line.strip()
    parts = line.split('_')
    try:
        #assume that everything that ends with '_number' is NCBI with GI
        gi = int(parts[-1])
        taxid = giToTaxa[gi]
        nodeTop = taxaTree.getNode(taxid).whichSupernode(topNodes)
        out.write("%s %s %s %s\n" % (line,nodeTop.id,nodeTop.name,taxaTree.getNode(taxid).lineageStr()))
    except:
        out.write("%s 0 NONE\n" % (line,))
    if iRec % 10000 == 0:
        print "Processed %s records" % (iRec+1,)
    iRec+=1
inp.close()
out.close()

