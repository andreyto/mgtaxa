"""Graphics for SOM model."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.Svm import *
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

class LabelMapper:
    def __init__(self,taxaTree):
        self.taxaTree = taxaTree
        self.topNodes = [ taxaTree.getNode(id) for id in micVirTaxids ]
        
    def label(self,id):
        (idPref,idSuf) = id.rsplit('_',1)
        if idPref == 'NCBIVM':
            idTopNode = self.taxaTree.getNode(int(idSuf)).whichSupernode(self.topNodes).id
            if idTopNode == 2157:
                idTopNode = 2
            lab = "%s_%s" % (idPref,idTopNode)
        elif idPref.startswith('GSIOVIR'):
            lab = 'GSIOVIR'
        elif idPref == 'GSIOMIC_10239':
            lab = 'GSIOVIR'
        elif idPref == 'GSIOMIC_2157':
            lab = 'NCBIVM_2'
        else:
            lab = idPref
        return lab

    def color(self,lab):
        if lab == 'NCBIVM_2':
            col = 'red'
        elif lab == 'NCBIVM_10239':
            col = 'blue'
        elif lab == 'GSIOVIR':
            col = 'green'
        elif lab == 'GSIOMIC_2':
            col = 'cyan'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        return col

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
col = [ labMapper.color(l) for l in lab ]
print zip(lab,col)
#col = plib.colors.colorConverter.to_rgba_list(col)
color = dict(zip(lab,col))
colorMap = {}
for (id,lab) in idToLab.items():
    colorMap[id] = color[lab]

out = open(opt.outColors,'w')
#out.write("id,red,green,blue,alpha\n")
out.write("id,color\n")
for id in ids:
    c = colorMap[id]
    #out.write("%s,%s,%s,%s\n" % (id,c[0],c[1],c[2],c[3]))
    out.write("%s,%s\n" % (id,c))
out.close()

