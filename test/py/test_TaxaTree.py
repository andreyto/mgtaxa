### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.TaxaTree import *
from MGT.TaxaIO import *

from MGT.Config import options

from time import time

def test_merged(tree):
    merged = tree.getMerged()
    id_merged = merged.keys()[0]
    node = tree.getNode(id_merged)
    assert merged[id_merged] == node.id

store = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,
        ncbiNamesDumpFile=options.taxaNamesFile,
        ncbiMergedDumpFile=options.taxaMergedFile)

print "Creating the tree from NCBI dump file"

start = time()

if False:
    def debugOnNodeUpdate(node,name,value):
        if getattr(node,'id',0) == 144409:
            print node.id, name, value

    TaxaTree.setDebugOnUpdate(debugOnNodeUpdate)
    taxaTree = TaxaTree(storage=store)
    taxaTree.unsetDebugOnUpdate()

    print "DEBUG: Done in %s sec" % (time() - start)

#for node in taxaTree.iterDepthTop():
#    print node

    taxaTree.deleteNodesIf(lambda n: n.rank == "species")
    taxaTree.reindex()
    taxaTree.setIndex()

else:
    taxaTree = TaxaTree(storage=store)
    
    print "DEBUG: Done in %s sec" % (time() - start)

test_merged(taxaTree)


#
#storeOut = NodeStorageNewick(fileName="test_TaxaTree.newick",
#    labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))

#storeOut = NodeStorageHypView(fileName="test_TaxaTree.hv3",
#    labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))

#storeOut = NodeStorageLibSea(fileName="test_TaxaTree.libsea",
#   labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))
        
#storeOut.save(taxaTree)


######### JSON #########

storeJson = NodeStorageJson(fileName="test_TaxaTree.json")


print "DEBUG: JSON dumping the tree"

start = time()

storeJson.save(taxaTree)

print "DEBUG: Done in %s sec" % (time() - start)


print "DEBUG: Creating the tree from JSON storage"

start = time()

taxaTreeJson = TaxaTree(storage=storeJson)

print "DEBUG: Done in %s sec" % (time() - start)

taxaTreeJson.reindex()

test_merged(taxaTreeJson)

####### Pickle #########

storePickle = NodeStoragePickle(fileName="test_TaxaTree.pkl")

print "DEBUG: Pickling the tree"

start = time()

storePickle.save(taxaTree)

print "DEBUG: Done in %s sec" % (time() - start)

taxaTree = None

print "DEBUG: Creating the tree from pickled storage"

start = time()

taxaTreePickle = TaxaTree(storage=storePickle)

print "DEBUG: Done in %s sec" % (time() - start)

taxaTreePickle.reindex()

#for node in taxaTreePickle.getNodesIter():
#    print node

