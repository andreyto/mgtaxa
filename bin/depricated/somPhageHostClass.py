### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Graphics for SOM model."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.GHSOM import *
from MGT.SOMCode import *
from MGT.SOMGraph import *
from MGT.Graphics import *
from MGT.Svm import *
from MGT.Phage import *
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-name",
        action="store", type="string",dest="inName"),
        make_option("-s", "--in-samp",
        action="append", type="string",dest="inSamp"),
        make_option("-p", "--in-ph-picks",
        action="store", type="string",dest="inPhPicks"),
        make_option("-d", "--mod-dump",
        action="store", type="string",dest="modDump"),
        make_option("-n", "--test-null",
        action="store_true", dest="testNull",default=False),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


def getNCBINode(taxaTree,idSamp):
    return taxaTree.getNode(int(idSamp))


opt,args = getProgOptions()

taxaTree = loadTaxaTreeNew()

micNodes = [ taxaTree.getNode(id) for id in micTaxids ]
#phageNode = taxaTree.getNode(phageTailedTaxid)
virNode = taxaTree.getNode(virTaxid)

if not os.path.isfile(opt.modDump):
    som = GHSOM(opt.inName) #SOMCode(opt.inName)
    mod = som.loadModel()
    samp = n.concatenate([ loadSparseSeqsAsDense(inSamp) for inSamp in opt.inSamp ])
    mod.setSamples(samp)
    print "Mapping samples..."
    #mod.mapSamples(pmRadius=mod.paretoRadiusLenNormSamp)
    # Pareto radius is too large - U*-matrix is dominated by P-matrix
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
unit = mod.unit
sampNode = mod.sampNode
samp = mod.getSamples()
nbCnt = {}
nodeCnt = {}

if opt.testNull:
    nodes = [ (id,getNCBINode(taxaTree,id)) for id in samp["id"] ]
    ids = set([ id for (id,node) in nodes if node is not None and node.isSubnode(virNode) ])
    ids = n.asarray(list(ids),dtype='O')
    idsRnd = permuteObjArray(ids)
    idRndMap = dict(zip(ids,idsRnd))
    for iSamp in xrange(len(samp)):
        id = samp[iSamp]["id"]
        if id in idRndMap:
            samp[iSamp]["id"] = idRndMap[id]

for iSamp in xrange(len(samp)):
    sm = samp[iSamp]
    id = sm["id"]
    #label = idToLab[samp[iSamp]["id"]]
    #if label == "NCBI Myoviridae"
    node = getNCBINode(taxaTree,id)
    #if node is not None and node.whichSupernode(micNodes) is not None:
    #    print "Mic node: ", node.lineageStr()
    try:
        nodeCnt[node] += 1
    except KeyError:
        nodeCnt[node] = 1
    if node is not None and node.isSubnode(virNode):
        cellInd = sampNode[iSamp]
        for nbSize in (0,):
            nbCnt1 = {}
            cells,cellsInd = getRectStencil(unit,cellInd,nbSize) #5
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
                            nbCnt1[pair] += 1
                        except KeyError:
                            nbCnt1[pair]  = 1
            for pair in nbCnt1:
                try:
                    nbCnt[pair] += 1
                except KeyError:
                    nbCnt[pair] = 1
            if len(nbCnt1) > 0:
                break
aCnt = {}
for (pair,cnt) in nbCnt.items():
    try:
        aCnt[pair[0]].append((cnt,pair[1]))
    except KeyError:
        aCnt[pair[0]] = [ (cnt,pair[1]) ]
for val in aCnt.values():
    val.sort(reverse=True)

for node,matches in aCnt.items():
    if float(matches[0][0])/nodeCnt[node] < 0:
        del aCnt[node]

phost = PhageHostSeqPicker(taxaTree=taxaTree)
phost.load(opt.inPhPicks)
virHosts = phost.seqVirHosts()

sampNodes = set([getNCBINode(taxaTree,id) for id in set(samp["id"])])
cntTP = 0
cntP = 0
print "Host assignments for each virus\n"
for (nodeName,node,matches) in sorted([ (node.name,node,matches) for (node,matches) in aCnt.items()]):
    nodeHosts = set(virHosts[node])
    #predHosts = [ (nodeM.name,cnt,nodeM in nodeHosts) for (cnt,nodeM) in matches ]
    predHost = matches[0][1] 
    cntTP +=  predHost in nodeHosts #or sum([h.lcsNode(predHost).rank == "family" for h in nodeHosts])>0
    #count this as prediction only if there are at least some host samples in the dataset
    hasHostSamples = len(nodeHosts & sampNodes) > 0
    if hasHostSamples:
        cntP += 1
    print "Phage: %s, host samples: %s" % (node.name,hasHostSamples)
    print "    Microbes: %s\n" % ", ".join(["%s : %s : %s" % (nodeM.name,cnt,nodeM in nodeHosts) for (cnt,nodeM) in matches])
print "TP: %s  P: %s SPE: %.0f TH: %s SPE_TH: %.0f" % (cntTP,len(aCnt),float(cntTP)/len(aCnt)*100,cntP,float(cntTP)/cntP*100)

