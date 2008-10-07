"""Graphics for SOM model."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.GHSOM import *
from MGT.SOMGraph import *
from MGT.Graphics import *
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

