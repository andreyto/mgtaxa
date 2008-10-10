"""Graphics for SOM model."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.GHSOM import *
from MGT.SOMGraph import *
from MGT.Graphics import *
from MGT.Svm import *
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-name",
        action="store", type="string",dest="inName"),
        make_option("-s", "--in-samp",
        action="store", type="string",dest="inSamp"),
        make_option("-d", "--mod-dump",
        action="store", type="string",dest="modDump"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

mgtDbDir = "/home/atovtchi/work/mgtdata"

taxaDir=os.path.join(mgtDbDir,"taxonomy.new")

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))

if not os.path.isfile(opt.modDump):
    som = GHSOM(opt.inName)
    som.setModelDir('.')
    mod = som.loadModel(components=("unit","weights"))
    samp = loadSparseSeqsAsDense(opt.inSamp)
    mod.setSamples(samp)
    print "Mapping samples..."
    #mod.mapSamples(pmRadius=mod.paretoRadiusLenNormSamp)
    mod.mapSamples(pmRadius=0.2) #0.2
    #mod.mapSamples(pmRadius=None)
    mod.makeData()
    dumpObj(mod,opt.modDump)
else:
    mod = loadObj(opt.modDump)
mod.makeUStarMatrix()
mod.makeUnit()
#umat = mod.umat
#mod.makeUMatrixTest()
#umat_test = mod.umat
#assert n.allclose(umat,umat_test)
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
colorConv = plib.colors.colorConverter
try:
    colToRgba = colorConv.to_rgba_list
except AttributeError:
    colToRgba = colorConv.to_rgba_array
col = colToRgba(col)
color = dict(zip(lab,col))
#colorMap = {}
#for (id,lab) in idToLab.items():
#    colorMap[id] = color[lab]
mod.setLabels(idToLab)
pl.gray()
plot = SOMPlot(mod=mod,labColorMap=color,markerSize=40,figSizeScatter=12,figSizePie=5,figSizeMat=12)
plot.plot()


