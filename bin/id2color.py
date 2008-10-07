"""Graphics for SOM model."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.Svm import *
from MGT.Graphics import *
#import pylab as pl
#import matplotlib as plib
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-ids",
        action="store", type="string",dest="inIds"),
        make_option("-o", "--out-colors",
        action="store", type="string",dest="outColors"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

mgtDbDir = "/home/atovtchi/work/mgtdata"

taxaDir=os.path.join(mgtDbDir,"taxonomy.new")

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))

ids = svmLoadId(opt.inIds)
idToLab = {}
labMapper = LabelMapper(taxaTree=taxaTree)
for id in ids:
    idToLab[id] = labMapper.label(id)
lab = sorted(set(idToLab.values()))
countsLab = binCount(idToLab.values())
col = [ labMapper.color(l) for l in lab ]
cnt = [ countsLab[l] for l in lab ]
cntSum = sum(countsLab.values())
cntPer = [ x*100./cntSum for x in cnt ]
print zip(lab,col,cnt,cntPer)
#col = plib.colors.colorConverter.to_rgba_list(col)
color = dict(zip(lab,col))
colorMap = {}
for (id,lab) in idToLab.items():
    colorMap[id] = color[lab]
pdb.set_trace()
out = open(opt.outColors,'w')
#out.write("id,red,green,blue,alpha\n")
out.write("id,color\n")
for id in ids:
    c = colorMap[id]
    #out.write("%s,%s,%s,%s\n" % (id,c[0],c[1],c[2],c[3]))
    out.write("%s,%s\n" % (id,c))
out.close()

