from MGT.Taxa import *
from MGT.Debug import *

t = Timer()

def collectLeafIds(topNode,dstAttr): 
    """assigns to dstAttr a list of all leaf node IDs under each node"""
    topNode.setReduction(lambda node: [ node.id ] if node.isLeaf() else [],dstAttr,operator.add)

def collectLeafIdsBounded(topNode,dstAttr,nMaxPicked):
    """same as collectLefIds(), but propagates only a max 'n' elements from a list accumulated for each child node
    (can be viewed as a bounded sampling of the tree)"""
    boundedSample = lambda l,n: random.sample(l,min(len(l),n))
    topNode.setReduction(lambda node: [ node.id ] if node.isLeaf() else [],dstAttr,operator.add,
            childExtractor=lambda child: boundedSample(getattr(child,dstAttr),nMaxPicked))

def collectLeafIdsBalanced(topNode,dstAttr,nMaxPicked):
    """same as collectLefIds(), but propagates only a max 'n' elements from a list accumulated for each child node
    (can be viewed as a bounded sampling of the tree)"""
    boundedSample = lambda l,n: random.sample(l,min(len(l),n))
    topNode.setReduction(lambda node: [ node.id ] if node.isLeaf() else [],dstAttr,operator.add,
            childExtractor=lambda child: boundedSample(getattr(child,dstAttr),nMaxPicked))

def testReduction(taxaTree):
    rootNode = taxaTree.getRootNode()
    nLeaf = 0
    for node in rootNode.iterDepthTop():
        if node.isLeaf():
            nLeaf += 1
    print t.msg("Leaves counted")    
    print Memuse.resident()
    collectLeafIds(rootNode,"tmp_leafIds")
    print t.msg("collectLeafIds")
    print Memuse.resident()
    assert len(rootNode.tmp_leafIds) == nLeaf

taxaTree = loadTaxaTree()
print t.msg("Tree loaded")
#This always gives 0 on CentOS 5:
#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

print Memuse.resident()
testReduction(taxaTree)

