from MGT.Taxa import *

taxaTree = loadTaxaTreeTest()

taxaLevels = TaxaLevels()
taxaLevels.setLevels(taxaTree)

def _printFixedLin(node,lin):
    print "------------------------------------------------"
    print "Actual  lineage:     %s" % (node.lineageStr(),)
    print
    print "Linnean lineage:     %s" % \
            (' <<->> '.join([ "rank=%s : name=\"\"\"%s\"\"\" : taxid=%s" % \
            (rank,nd.name if nd else None,nd.id if nd else None) \
            for (rank,nd) in zip(taxaLevels.getLevelNames("ascend"),lin)]),)
    print "------------------------------------------------"

def _getFixedLinAndPrint(node):
    _printFixedLin(node,taxaLevels.lineageFixedList(node,null=None,format="node",fill="up-down"))


node = taxaTree.getNode(62672) #SAR86 cluster - no_rank, species under; no_rank, then class above
_getFixedLinAndPrint(node)
node = taxaTree.getNode(159775) #uncultured alpha proteobacterium MB11B03, species, several no_ranks above
_getFixedLinAndPrint(node)
node = taxaTree.getNode(62654) #SAR 116 cluster - no_rank, genus and species under; no_rank, then class above
_getFixedLinAndPrint(node)
node = taxaTree.getNode(883078) #Afipia broomeae ATCC 49717 - no_rank (strain), all levels present
_getFixedLinAndPrint(node)

#taxaLevels.reduceNodes(taxaTree)

#for node in taxaTree.getNodesIter():
#    print node
