### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Time the creation of many python objects that are used as tree nodes."""

from MGT.Taxa import *

from MGT.Config import options

store = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=None)

#nodes = store.load()

#print nodes[1]

taxaTree = TaxaTree(storage=store)
print taxaTree.getNode(1239)
