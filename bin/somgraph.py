"""Graphics for SOM model."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.GHSOM import *
from MGT.SOMGraph import *
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-name",
        action="store", type="string",dest="inName"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args

class LabelMapper:
    def __init__(self,taxaTree):
        self.taxaTree = taxaTree
        self.topNodes = [ taxaTree.getNode(id) for id in (phageTailedTaxid,)+micVirTaxids ]
        
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
            col = 'r'
        elif lab == 'NCBIVM_10239':
            col = 'b'
        elif lab == 'NCBIVM_28883':
            col = 'c'
        elif lab == 'GSIOVIR':
            col = 'g'
        elif lab == 'GSIOMIC_2':
            col = 'y'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        return col

opt,args = getProgOptions()

mgtDbDir = "/home/andrey/work/mgtdata"

taxaDir=os.path.join(mgtDbDir,"taxonomy.new")

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))

som = GHSOM(opt.inName)
som.setModelDir('.')
mod = som.loadModel()
idsSom = mod.sampIds()
idToLab = {}
labMapper = LabelMapper(taxaTree=taxaTree)
for id in idsSom:
    idToLab[id] = labMapper.label(id)
lab = sorted(set(idToLab.values()))
#print lab
#col = n.arange(len(lab),dtype=float) + 10
col = [ labMapper.color(l) for l in lab ]
print zip(lab,col)
#pdb.set_trace()
#cmap = pl.cm.spectral
#cnorm = pl.normalize(0,col.max())
#col = cmap(cnorm(col))
col = plib.colors.colorConverter.to_rgba_list(col)
color = dict(zip(lab,col))
#colorMap = {}
#for (id,lab) in idToLab.items():
#    colorMap[id] = color[lab]
mod.setLabels(idToLab)
somPlotScatter(mod,colorMap=color,markerSize=40,figSizeScatter=(12,12),figSizePie=(8,8))

