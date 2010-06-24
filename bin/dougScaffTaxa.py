### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Load and pickle Doug's scaffold-level taxonomy assignment file."""
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

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))

nameToId = taxaTree.makeNameToIdMap()

nameToIdTr = dict( ( (key.replace(" ","_"),val) for (key,val) in nameToId.items() ) )

inp = openCompressed(opt.inFile,'r')

fldDtype=[("id","i8"),("length","i8"),("level","i4"),("num_pep","i4"),("taxid","i4"),("per_pep","f4")]
fldValues = []

for row in inp:
    parts = row.split()
    name = parts[-2]
    try:
        taxid = nameToIdTr[name]
        parts[-2] = taxid
        parts = [ float(p) for p in parts ]
        fldValues.append(tuple(parts))
    except KeyError:
        #print "Name %s not found in taxonomy database, record skipped: %s" % (name,row)
        pass

res = numpy.asarray(fldValues,dtype=fldDtype)
dumpObj(res,opt.outFile)

