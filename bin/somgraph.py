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
topNodes = [ taxaTree.getNode(id) for id in micVirTaxids ]
for id in idsSom:
    (idPref,idSuf) = id.rsplit('_',1)
    if idPref == 'NCBIVM':
        idTopNode = taxaTree.getNode(int(idSuf)).whichSupernode(topNodes).id
        lab = "%s_%s" % (idPref,idTopNode)
    elif idPref.startswith('GSIOVIR'):
        lab = 'GSIOVIR'
    else:
        lab = idPref
    idToLab[id] = lab
lab = sorted(set(idToLab.values()))
color = dict(zip(lab,n.arange(len(lab))))
colorMap = {}
for (id,lab) in idToLab.items():
    colorMap[id] = color[lab]
mod.setLabels(colorMap)
somPlotScatter(mod,markerSize=16)

