from MGT.Taxa import *

from MGT.Config import options

taxaTree = TaxaTree(ncbiDumpFile=options.taxaNodesFile,save=False,load=False)

for node in taxaTree.iterDepthTop():
    print node

taxaLevels = TaxaLevels(taxaTree)

for node in taxaTree.getNodesIter():
    print taxaLevels.lineageKeys(node), node
