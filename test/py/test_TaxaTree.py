from MGT.TaxaTree import *
from MGT.TaxaIO import *

from MGT.Config import options

store = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile)

taxaTree = TaxaTree(storage=store)

#for node in taxaTree.iterDepthTop():
#    print node

taxaTree.deleteNodesIf(lambda n: n.rank == "species")
taxaTree.reindex()


#taxaLevels = TaxaLevels()
#taxaLevels.setLevels(taxaTree)
#taxaLevels.reduceNodes(taxaTree)

#for node in taxaTree.getNodesIter():
#    print node

#storeOut = NodeStorageNewick(fileName="test_TaxaTree.newick",
#    labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))

#storeOut = NodeStorageHypView(fileName="test_TaxaTree.hv3",
#    labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))

storeOut = NodeStorageLibSea(fileName="test_TaxaTree.libsea",
   labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))
        
storeOut.save(taxaTree)
