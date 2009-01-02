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
        action="append", type="string",dest="inSamp"),
        make_option("-d", "--mod-dump",
        action="store", type="string",dest="modDump"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


def getNCBINode(taxaTree,idSamp):
    (idPref,idSuf) = idSamp.rsplit('_',1)
    if idPref.startswith('NCBI'):
        return taxaTree.getNode(int(idSuf))
    else:
        return None


opt,args = getProgOptions()

taxaTree = loadTaxaTreeNew()

micNodes = [ taxaTree.getNode(id) for id in micTaxids ]
phageNode = taxaTree.getNode(phageTailedTaxid)

if not os.path.isfile(opt.modDump):
    som = GHSOM(opt.inName)
    som.setModelDir('.')
    mod = som.loadModel(components=("unit","weights"))
    samp = n.concatenate([ loadSparseSeqsAsDense(inSamp) for inSamp in opt.inSamp ])
    mod.setSamples(samp)
    print "Mapping samples..."
    #mod.mapSamples(pmRadius=mod.paretoRadiusLenNormSamp)
    # Pareto radius is too large - U*-matrix is domnated by P-matrix
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
### print taxa of samples close to myoviridae
unit = mod.unit
sampNode = mod.sampNode
samp = mod.getSamples()
nbCnt = {}
for iSamp in xrange(len(samp)):
    sm = samp[iSamp]
    id = sm["id"]
    #label = idToLab[samp[iSamp]["id"]]
    #if label == "NCBI Myoviridae"
    node = getNCBINode(taxaTree,id)
    #if node is not None and node.whichSupernode(micNodes) is not None:
    #    print "Mic node: ", node.lineageStr()
    if node is not None and node.isSubnode(phageNode) and not node.name.startswith("Enterobacteria"):
        cellInd = sampNode[iSamp]
        cells,cellsInd = getRectStencil(unit,cellInd,5)
        for cell in cells.flat:
            for idNb in cell:
                nodeNb = getNCBINode(taxaTree,idNb)
                if nodeNb is not None and nodeNb.whichSupernode(micNodes) is not None:
                    print "Phage: ", node.name #node.lineageStr()
                    #print "\n"
                    print "Microbe: ", nodeNb.name #nodeNb.lineageStr()
                    print "\n\n\n"
                    pair = (node,nodeNb)
                    try:
                        nbCnt[pair] += 1
                    except KeyError:
                        nbCnt[pair]  = 1
aCnt = {}
for (pair,cnt) in nbCnt.items():
    try:
        aCnt[pair[0]].append((cnt,pair[1]))
    except KeyError:
        aCnt[pair[0]] = [ (cnt,pair[1]) ]
for val in aCnt.values():
    val.sort(reverse=True)
print "Match counts for each phage\n"
for (node,matches) in aCnt.items():
    print "Phage: ", node.name
    print "    Microbes: %s\n" % ", ".join(["%s:%s" % (nodeM.name,cnt) for (cnt,nodeM) in matches])
sys.exit(0)
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
def onClick(event):
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%\
            (event.button, event.x, event.y, event.xdata, event.ydata)
    ids = mod.unit[round(event.xdata),round(event.ydata)]
    print ids
    for id in ids:
        if id.startswith('NCBIVM_'):
            taxid = int(id.split('NCBIVM_')[1])
            node = taxaTree.getNode(taxid)
            print node.lineageStr()
pl.gray()
plot = SOMPlot(mod=mod,labColorMap=color,markerSize=40,figSizeScatter=12,figSizePie=5,figSizeMat=12,onClick=onClick)
plot.plot()
pl.show()

