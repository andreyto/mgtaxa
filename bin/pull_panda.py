### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Pull sequence for a given set of taxa directly from PANDA fasta file.
Created to test the correctness of MGTAXA own sequence/sample database.
The current implementation collects everything in memory."""

from MGT.Common import *
from itertools import izip
from MGT.Taxa import *


pandaVirF = "/usr/local/db/panda/NuclSequences/ViralGroupNucl.fasta"

giToTaxa = loadObj("/home/atovtchi/work/mgtdata/taxa.pkl.gz")

#newTaxaDir="/home/atovtchi/work/mgtdata/taxonomy"

#taxaTreeNew = loadTaxaTree(ncbiDumpFile=os.path.join(newTaxaDir,"nodes.dmp"),
#        ncbiNamesDumpFile=os.path.join(newTaxaDir,"names.dmp"))
#taxaTree = loadTaxaTree()

inpSeq = FastaReader(pandaVirF)

iRec = 0
taxMis = Struct(bounds=0,zeroP=0,zeroN=0,nonzero=0)

for rec in inpSeq.records():
    hdr = rec.header()
    try:
        taxid = int(hdr.split("taxon:")[1].split(None,1)[0])
    except (ValueError,IndexError):
        taxid = 0
    gi = int(hdr.split("|",2)[2].split(None,1)[0])
    taxMisFlag = False
    if len(giToTaxa) <= gi:
        taxMis.bounds += 1
        taxMisFlag = True
    else:
        taxidGi = giToTaxa[gi]
        if taxidGi != taxid:
            if taxidGi == 0:
                taxMis.zeroN += 1
            elif taxid == 0:
                taxMis.zeroP += 1
            else:
                taxMis.nonzero += 1
            taxMisFlag = True
    if taxMisFlag:
        print taxid,hdr
    iRec += 1

print "Scanned %d records, taxonomy mismatches are:\n%s" % (iRec,taxMis)
inpSeq.close()

