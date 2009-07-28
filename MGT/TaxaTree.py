"""Classes to represent and traverse a taxonomy tree (NCBI format)"""

from MGT.Common import *
from MGT.TaxaConst import *

import os
import numpy
from numpy import array
import itertools, operator
from copy import copy
        
class TaxaNode(object):
    """Node of a taxonomy tree."""

    # This reduces memory consumption by removing __weakref__ dict attribute
    __slots__ = '__dict__'

    @staticmethod
    def setDebugOnUpdate(hook):
        """Sets new debugging hook to call on attribute update.
        @param hook - unbound function with signature f(node,name,value) where
        name is attribute name and value is its new value.
        @return previous hook function (or None), which can be stored by the caller and
        used later to restore the original behavior."""
        if hasattr(TaxaNode,'_hookDebugOnUpdate'):
            orig = TaxaNode._hookDebugOnUpdate
        else:
            orig = None
        if hook is None:
            TaxaNode.unsetDebugOnUpdate()
        else:
            TaxaNode._hookDebugOnUpdate = hook
            TaxaNode.__setattr__ = TaxaNode._debugSetAttr
        return orig

    @staticmethod
    def unsetDebugOnUpdate():
        """Remove checks for debugging on attribute update previously set by setDebugOnUpdate().
        This removes any run-time penalty for the debugging hook."""
        if hasattr(TaxaNode,'_hookDebugOnUpdate'):
            delattr(TaxaNode,'_hookDebugOnUpdate')
            delattr(TaxaNode,'__setattr__')

    def __init__(self):
        self.children = []
    
    def __str__(self):
        return strAttributes(self,exclude=("children","par"))

    def __getstate__(self):
        """Pickling support - exclude node cross references to improve memory and performance 100x.
        The TaxaTree constructor can restore the internal tree structure from remaining data.
        The TaxaTreeIO.NodeStoragePickle relies on this."""
        # The copy/del combination is a little faster than the iter comprehension
        #return dict((item for item in self.__dict__.iteritems() if item[0] not in ('par','children')))
        # Shallow copy
        newdict = copy(self.__dict__)
        del (newdict['par'], newdict['children'])
        return newdict

    def _debugSetAttr(self,name,value):
        """This method will replace the default __setattr__ and call the debugging hook.
        This is set by by @see setDebugOnUpdate()."""
        self._hookDebugOnUpdate(name,value)
        object.__setattr__(self,name,value)

    def setParent(self,par):
        """Double-link this node and the parent node.
        For speed, we do not check that 'par' is not a parent of this node already,
        so this method should not be called twice."""
        #TODO: is garbage collector going to deal with circular refs?
        self.par = par
        if par is not None:
            par.children.append(self)
        
    def removeChild(self,child):
        self.children.remove(child)
        
    def getParent(self):
        return self.par
    
    def getParentId(self):
        return self.idpar
    
    def getChildren(self):
        return self.children
        
    def numChildren(self):
        return len(self.children)

    def isRoot(self):
        return self.par is None
    
    def isLeaf(self):
        return len(self.children) == 0
    
    def isInternal(self):
        return not (self.par is None or len(self.children) == 0)

    def lineage(self,topNode=None):
        """Return a list of nodes comprising this nodes' lineage (from bottom to top, including this node).
        Setting the value of 'topNode' to some existing node allows to stop the lineage construction
        at the root of some sub-tree. If some node that is not in the lineage of this node is 
        provided as 'topNode', the function will raise an AttributeError when it reaches the root node
        and tries to access its parent."""
        node = self
        lin = [node]
#        linset = set()
        while node is not topNode:
# Debugging code
#            if node in linset:
#                print "Circular tree reference in 'lineage()': self = %s\n node = %s\n, lin = %s" % (self.data,node.data,(n.id for n in lin))
#                raise ValueError()
#            linset.add(node)
            node = node.par
            if node is not None:
                lin.append(node)
        return lin

    def lineageWhile(self,condition):
        """Return lineage that is cut as soon as an attempt is made to look above the root note or condition(node) returns False.
        Can possibly return an empty list.
        Condition is always called with a valid node argument."""
        node = self
        lin = []
        while node is not None and condition(node):
            lin.append(node)
            node = node.par
        return lin

    def findRankInLineage(self,rank,topNode=None):
        node = self
        if node.rank == rank:
            return node
        while node is not topNode:
            node = node.par
            if node is not None:
                if node.rank == rank:
                    return node
        return None
    
    def findLeftRankInLineage(self,ranks,topNode=None):
        """Return a node corresponding to the left-most element from 'ranks' found in the lineage"""
        ranksMap = dict([ (n.rank,n) for n in self.lineage(topNode=topNode) ])
        for rank in ranks:
            if rank in ranksMap:
                return ranksMap[rank]
        return None

    def lineageRanks(self,*l,**kw):
        return [ node.rank for node in self.lineage(*l,**kw) ]
    
    def lineageRanksTaxa(self,*l,**kw):
        return [ (node.rank,node.id) for node in self.lineage(*l,**kw) ]
    
    def lineageRanksStr(self,*l,**kw):
        return ','.join([ "%s=%s" % (rank,taxid) for (rank,taxid) in self.lineageRanksTaxa(*l,**kw) ])

    def lineageStr(self,*l,**kw):
        return ' <<->> '.join([ "rank=%s : name=\"\"\"%s\"\"\" : taxid=%s" % (node.rank,node.name,node.id) for node in self.lineage(*l,**kw) ])

    def visitDepthTop(self,func):
        """Apply function 'func' to each node traversing the subtree depth-first, applying 'func' before visiting the children.
        If function 'func' returns anything when applied to this node, stop the iteration and return w/o proceeding to the children."""
        if func(self) is None:
            for child in self.children:
                child.visitDepthTop(func)
    
    def iterDepthTop(self):
        """Return depth-first top-to-bottom iterator."""
        yield self
        for child in self.children:
            for node in child.iterDepthTop():
                yield node


    def iterDepthTwice(self):
        """Return depth-first on-entry+on-exit iterator (modified preorder).
        Dereferencing the iterator yields (node,visit) tuple, where visit is either 1 or 2 according to the time that we visit this node.
        See setNestedSetsIndex() method for application example."""
        yield (self,1)
        for child in self.children:
            for ret in child.iterDepthTwice():
                yield ret
        yield (self,2)

    def iterFillUp(self,topNode=None,topFirst=False):
        """Iterate as though we were pouring water into this node until it fills a given super-tree above it.
        @param topNode - "fill up" up to this super-node (inclusive), None means up to root node.
        @param topFirst - if True, the top node of each super-tree will be iterated before its sub-nodes, otherwise - after ("self" is always iterated first)."""
        #assert self.isSubnode(topNode)
        for node in self.iterDepthTop():
            yield node
        lin = self.lineage(topNode=topNode)
        for sub,sup in zip(lin[:-1],lin[1:]):
            if topFirst:
                yield sup
            for supChild in sup.children:
                if supChild is not sub:
                    for node in supChild.iterDepthTop():
                        yield node
            if not topFirst:
                yield sup

    def visitDepthBottom(self,func):
        """Apply function 'func' to each node traversing the subtree depth-first, applying 'func' after visiting the children."""
        #print "DEBUG: enter: ", self, len(self.children)
        for child in self.children:
            child.visitDepthBottom(func)
        func(self)
        #print "DEBUG: exit: ", self

    def visitDepthTwice(self,func):
        """Apply function 'func' to each node twice traversing the subtree depth-first. This is also called modified pre-order traverse.
        'func' is applied first time when coming to the node, in which case it is called as func(node,visit=1),
        and second time before returning from the node, in which case  it is called as func(node,visit=2)."""
        func(self,visit=1)
        for child in self.children:
            child.visitDepthTwice(func)
        func(self,visit=2)


    def isUnclassified(self):
        #DEBUG:
        #return False
        return self.name.lower().startswith("unclassified")

    def getDepth(self):
        """Return distance to the top of the tree.
        Depth must be precalculated for each node by a call to TaxaTree.setDepth()."""
        return self.depth

    def setDepth(self,depth):
        """Recursively set depth for this node and its subtree."""
        self.depth = depth
        for child in self.children:
            child.setDepth(depth+1)

    def getNested(self):
        return (self.lnest,self.rnest)

    def setNestedSetsIndexByIter(self,startIndex=1):
        ind = startIndex
        for (node,visit) in self.iterDepthTwice():
            if visit == 1:
                node.lnest = ind
            else:
                node.rnest = ind
            ind += 1

    def _setNestedVar(self,cnt):
        self.lnest = cnt
        cnt += 1
        for child in self.children:
            cnt = child._setNested(cnt)
        self.rnest = cnt
        return cnt + 1

    def _setNested(self,cnt):
        self.lnest = cnt.next()
        for child in self.children:
            child._setNested(cnt)
        self.rnest = cnt.next()

    def setNestedSetsIndex(self,startIndex=1):
        self._setNested(itertools.count(startIndex))

    def setIndex(self,startIndex=0):
        """Index nodes in a subtree in a depth first traversal order.
        For each node in a subtree, set attribute lind is the index of the node itself, 
        attribute rind is the maximum plus one
        index of all nodes in a subtree, including the top node. Index starts at startIndex
        and only incremented when node is entered the first time in depth first traversal.
        Thus, for any subtree starting at some node, it is easy to create an attribute map as
        a dense (Numpy) array: a = numpy.zeros(top.rind - top.lind); a[node.lind-top.lind] = value.
        It also makes it easy and fast to work with any subtree below the node:
        sub_a = a[sub_top.lind-top.lind:sub_top.rind-top.lind]
        sub_a will densely list all values for nodes below the sub_top, with elements accessed as
        sub_a[node.lind-sub_top.lind]."""
        def actor(node,i):
            node.lind = i
            for child in node.children:
                i = actor(child,i+1)
            node.rind = i + 1
            return i
        actor(self,startIndex)

    def makePropArrayZeros(self,dtype):
        """Return (self.lind,numpy.zeros(number of nodes in subtree))."""
        return (self.lind,numpy.zeros(self.rind - self.lind,dtype=dtype))

    def setPropArrayFromAttrib(self,prop,name,base=None):
        """Assign values of subtree node attributes to numpy property array.
        @param prop - Numpy array
        @param name - node attribute name
        @param base - indexing will be done as prop[subnode.lind-base]. If 
        base is None, the lind of this node will be used."""
        if base is None:
            base = self.lind
        for node in self.iterDepthTop():
            prop[node.lind-base] = getattr(node,name)

    def setPropArrayFromFunc(self,prop,func,base=None):
        """Assign return value of a func(subnode) to numpy property array.
        @param prop - Numpy array
        @param func - unary function acting on the subtree node
        @param base - indexing will be done as prop[subnode.lind-base]. If 
        base is None, the lind of this node will be used."""
        if base is None:
            base = self.lind
        for node in self.iterDepthTop():
            prop[node.lind-base] = func(node)
    
    def setAttribFromPropArray(self,prop,name,base=None):
        if base is None:
            base = self.lind
        for node in self.iterDepthTop():
            # a[i].item() converts to Python type,
            # Numpy docs say it is supposed to be faster 
            # than Numpy scalar, and it indeed speeds up
            # tree operations by x10
            setattr(node,name,prop[node.lind-base].item())

    def makePropArrayFromAttrib(self,dtype,name):
        (base,prop) = self.makePropArrayZeros(dtype=dtype)
        self.setPropArrayFromAttrib(prop=prop,name=name,base=base)
        return prop

    def makePropArrayFromFunc(self,dtype,func):
        (base,prop) = self.makePropArrayZeros(dtype=dtype)
        self.setPropArrayFromFunc(prop=prop,func=func,base=base)
        return prop
    
    def setTotal(self,srcAttr,dstAttr):
        """Set total as 'dstAttr' attribute as a sum of 'srcAttr' attribute of this node and all its subnodes."""
        assert srcAttr != dstAttr
        
        def actor(node):
            setattr(node,dstAttr,sum((getattr(child,dstAttr) for child in node.children),getattr(node,srcAttr)))
        
        self.visitDepthBottom(actor)


    def setReduction(self,extractor,dstAttr,reduction=None,condition=None):
        """Assign to 'dstAttr' attribute the result of applying 'reduce(func,...)' to this node and its children.
        Example: setTotal() that sums only for non "unclassified" children could be implemented as:
        setReduction(lambda node: getattr(node,srcAttr),dstAttr,operator.add,lambda node: not node.isUnclassified()).
        Note: extractor(node) is called unconditionally; condition is only applied when the computed values of
        'dstAttr' attribute are accumulated for reduction.
        @param extractor - unary function object that is applied to each node to extract the value that is accumulated.
        As a convenient exception, if 'extractor' is a string, it is considered to be a node attribute name.
        @param dstAttr - name of destination attribute to set in each node with accumulated result
        @param reduction - binary function object that is passed as a first argument to 'reduce' built-in, defaults
        to operator.add
        @param condition - apply this condition to each child node to decide if it should contribute to the 'dstAttr'
        of the parent."""
        if isinstance(extractor,StringTypes):
            srcAttr = extractor
            extractor = lambda node: getattr(node,srcAttr)
        if reduction is None:
            reduction=operator.add
        if condition is None:
            condition=lambda node: True
        def actor(node):
            setattr(node,
                    dstAttr,
                    reduce(reduction,
                           (getattr(child,dstAttr) for child in node.children if condition(node)),
                           extractor(node)
                          )
                    )
        self.visitDepthBottom(actor)

    def setMaxSubtreeRank(self):
        """Set an attribute "rank_max" in every node of this subtree that is a highest standard rank of any node in a subtree.
        unclassRank and noRank nodes are both considered as noRank in comparisons.
        Thus, if and only if all subtree nodes are unclassRank or noRank, "rank_max" is set to noRank."""
        dstAttr = "rank_max"
        #{rank->order index}
        ranksToInd = dict([ (x[1],x[0]) for x in enumerate((noRank,)+tuple(linnRanks)) ])
        def extractor(node):
            """Return rank converted such that everything not found in ranksToInd becomes noRank.
            That includes unclassRank"""
            rank = node.rank
            if rank == unclassRank:
                return noRank
            elif rank in ranksToInd:
                return rank
            else:
                return noRank
        def reduction(rank_max,rank):
            if ranksToInd[rank_max] >= ranksToInd[rank]:
                return rank_max
            else:
                #because self.rank is supplied as starting value to reduce(),
                #it is weird if any standard rank for max turns out to be
                #lower than for a subnode
                assert self.rank in (noRank,unclassRank)
                return rank
        self.setReduction(extractor=extractor,dstAttr=dstAttr,reduction=reduction)

    def setAttribute(self,name,value):
        for node in self.iterDepthTop():
            setattr(node,name,value)

    def delAttribute(self,name):
        for node in self.iterDepthTop():
            try:
                delattr(node,name)
            except AttributeError:
                pass

    def isSubnode(self,other):
        """Return true if this node is a descendant of the other node.
        Uses pre-computed nested sets indexes, which must be up-to-date.
        A node is not considered a subnode of itself."""
        return self.lnest > other.lnest and self.rnest < other.rnest

    def isSubnodeAny(self,others):
        return self.whichSupernode(others) is not None

    def whichSupernode(self,others):
        """Return node from others that is a supernode of self, or None"""
        for o in others:
            if self.isSubnode(o):
                return o
        return None

    def lcsNode(self,other):
        """Return the Lowest Common SuperNode of self and other"""
        if self.getDepth() < other.getDepth():
            x = self
            y = other
        else:
            x = other
            y = self
        for node in x.lineage():
            if y.isSubnode(node) or y is node:
                return node
        assert ValueError, "Could not find lcsNode. Tree is corrupt or nodes are from differnt trees"

    def setIsUnderUnclass(self):
        """Set an attribute 'isUnderUnclass' for the nodes in this branch.
        It flags all nodes that have "unclassified" super-node somewhere in their lineage
        up to this node including. This node will be marked only if it isUnclassifed()."""
        iter = self.iterDepthTop()
        node = iter.next()
        node.isUnderUnclass = node.isUnclassified()
        for node in iter:
            if node.isUnclassified() or node.getParent().isUnderUnclass:
                node.isUnderUnclass = True
            else:
                node.isUnderUnclass = False


class TaxaTree(object):
    
    @staticmethod
    def setDebugOnUpdate(hook):
        """Calls @see TaxaNode.setDebugOnUpdate()"""
        return TaxaNode.setDebugOnUpdate(hook)


    @staticmethod
    def unsetDebugOnUpdate():
        """Calls @see TaxaNode.unsetDebugOnUpdate()"""
        TaxaNode.unsetDebugOnUpdate()

    def __init__(self,storage):
        """Take a storage instance (such as NodeStorageNcbiDump), extract node data from it, link nodes
        with python object references and setup other attributes if necessary.
        Node links are created based on getParentId() calls to each node. For the root node, and root node only,
        this method must return a value that evaluates to False in logical expression (e.g. 0 or None)."""

        nodes = storage.load()
        self.nodes = nodes
        self.rootNode = None
        for v in nodes.itervalues():
            #We use None for parent link of the root node
            try:
                par = nodes[v.idpar]
            except KeyError:
                par = None
                self.rootNode = v
            v.setParent(par)

        rootNode = self.getRootNode()
        assert rootNode is not None, "Root node was not detected in TaxaTree input"
        ## Calculate additional node attributes if the root node does not
        ## have them already
        try:
            rootNode.getDepth()
        except AttributeError:
            self.setDepth()
        try:
            rootNode.getNested()
        except AttributeError:
            self.setNestedSetsIndex()
    
    def write(self,out):
        for id in sorted(self.nodes.iterkeys()):
             out.write(str(self.nodes[id])+'\n')
    
    def writeLineage(self,out):
        for id in sorted(self.nodes.iterkeys()):
            out.write(self.nodes[id].lineageRanksStr()+'\n')
             
    def setDepth(self):
        self.getRootNode().setDepth(0)
            
    def setNestedSetsIndex(self,startIndex=1):
        self.getRootNode().setNestedSetsIndex(startIndex=startIndex)

    def reindex(self):
        """Refresh internal derived indices (depth, nested sets).
        You probably want to call it after changing the tree structure (e.g. deleting nodes).
        This method DOES NOT call setIndex()."""
        self.setDepth()
        self.setNestedSetsIndex()

    def setTotal(self,srcAttr,dstAttr):
        self.getRootNode().setTotal(srcAttr,dstAttr)

    def setReduction(self,extractor,dstAttr,reduction=None,condition=None):
        self.getRootNode().setReduction(extractor,dstAttr,reduction,condition)

    def getNode(self,id):
        return self.nodes[id]

    def getNodes(self,ids):
        return [ self.nodes[id] for id in ids ]

    def removeNodes(self,ids):
        """Remove nodes referenced by 'ids' from auxiliary data structures.
        @precondition - Internal parent node links should no longer reference these nodes.
        @postcondition - Indexing methods such as getNodesDict() and getNodesIter() will
        no longer list these nodes.
        This bookkeeping method should be called to maintain tree consistency after any re-arrangement is
        done to the tree structure."""
        nodes = self.nodes
        for id in ids:
            del nodes[id]
    
    def getNodesDict(self):
        return self.nodes
    
    def numNodes(self):
        return len(self.nodes)

    def numLinks(self):
        return self.numNodes() - 1
    
    def getNodesIter(self):
        return self.nodes.itervalues()
    
    def getRootNode(self):
        return self.rootNode

    def setAttribute(self,name,value,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        top.setAttribute(name,value)

    def delAttribute(self,name,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        top.delAttribute(name)

    def visitDepthTop(self,func,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        top.visitDepthTop(func)
        
    def iterDepthTop(self,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        return top.iterDepthTop()

    def iterDepthTwice(self,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        return top.iterDepthTwice()

    def visitDepthBottom(self,func,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        top.visitDepthBottom(func)

    def deleteNodesIf(self,predicate,id=None):
        """Delete every subtree which top node has 'predicate' == True"""
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        nodes = self.nodes
        for node in top.iterDepthTop():
            to_del = [child for child in node.getChildren() if predicate(child)]
            for child in to_del:
                node.removeChild(child)
                for n in child.iterDepthTop():
                    del nodes[n.id]

    def setIndex(self,startIndex=0,id=None):
        """See TaxaNode.setIndex() description."""
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        top.setIndex(startIndex=startIndex)


    def makeNameIndex(self):
        nameInd = {}
        namesInd = {}
        for node in self.iterDepthTop():
            name = node.name.lower()
            if name in nameInd:
                nameInd[name].append(node)
            else:
                nameInd[name] = [ node ]
            if hasattr(node,"names"):
                for name in [ n.lower() for n in node.names ]:
                    if name in namesInd:
                        namesInd[name].append(node)
                    else:
                        namesInd[name] = [ node ]
                    
        self.nameInd = nameInd
        self.namesInd = namesInd
        self.emptyList = []

    def searchName(self,name):
        try:
            return self.nameInd[name.lower()]
        except KeyError:
            return self.emptyList

    def searchNames(self,name):
        try:
            return self.namesInd[name.lower()]
        except KeyError:
            return self.emptyList

    def buildSubtreeMask(self,ids,values=None,other=-1):
        """Return a mask that 'paints' every node with a value from 'values' according to if a node is in a subtree
        of corresponding id from 'ids'. If 'values' is None, values of 'ids' are used.
        Mask for a non-existing node id is set to 0; Mask for nodes not from selected subtrees is set to 'other' value.
        'ids' will be processed in the order of elements (important if they point to overlapping subtrees)."""
        mask = self.presenceMask(value=other)
        if values is None:
            values = ids
        class mask_setter:
            def __init__(self,val):
                self.val = val
            def __call__(self,node):
                mask[node.id] = self.val
        for (id,val) in itertools.izip(ids,values):
            self.visitDepthTop(mask_setter(val=val),id=id)
        return mask

    def buildRankSubtreeMask(self,rank,other=-1):
        """Return a mask that 'paints' every node with an id of its closest node above that has rank 'rank'.
        Mask for a non-existing node id is set to 0; Mask for nodes not from selected subtrees is set to 'other' value."""
        mask = self.presenceMask(value=other)
        def mask_setter(node):
            node_id = node.id
            node_rank = node.rank
            if node_rank == rank:
                mask[node_id] = node_id
            elif node.par is not None:
                mask[node_id] = mask[node.par.id]
        self.visitDepthTop(mask_setter)
        return mask

    def makeNameToIdMap(self):
        x = {}
        for node in self.iterDepthTop():
            x[node.name] = node.id
        return x

    def maxId(self):
        try:
            return self._max_id
        except AttributeError:
            self.buildIdPresenceMask()
            return self._max_id

    def buildIdPresenceMask(self):
        ids = numpy.array(self.nodes.keys(),dtype=int)
        self._max_id = numpy.max(ids)
        mask = numpy.zeros(self.maxId()+1,dtype=bool)
        mask[ids] = 1
        self._mask_presence = mask
        
    def emptyMask(self,dtype=int):
        return numpy.zeros(self.maxId()+1,dtype=dtype)
    
    def presenceMask(self,value=1,dtype=int):
        try:
            mask = self._mask_presence
        except AttributeError:
            self.buildIdPresenceMask()
            mask = self._mask_presence
        return (mask * value).astype(dtype)


def selectTopNodes(nodes):
    """Select nodes that are not descendants of any other nodes in the input list.
    @param sequence of TaxaNodes
    @ret maximum subset of nodes such as not x.isSubnode(y) for any x,y.
    Average complexity is currently quadratic in len(nodes) (stupid implementation)
    Worst case complexity is always qudratic (e.g. if this method is called on its own result)
    """
    isSub = [False]*len(nodes)
    for i in xrange(1,len(nodes)):
        x = nodes[i]
        for j in xrange(0,i):
            y = nodes[j]
            if x.isSubnode(y):
                isSub[i] = True
            elif y.isSubnode(x):
                isSub[j] = True
    return [ nodes[i] for i in xrange(len(nodes)) if not isSub[i] ]


"""
total length of taxa
agregate original tree along (phylum, class, family etc) levels -
    leaves are the same (original taxa, except maybe strains are 
    agregated up to species). call AT - agregated tree, OT - original tree.
    For viruses and phages
bottom up assign to each node (AT or OT) the sum of all node lengths below it (leaves'
    lengths are their total sequence length expressed in the number of k-mer vectors,
    e.g. number of chunks of 50kb if we calculated vectors on 50kb chunks).
The goal: we need to asign the amount of sequence selected for a given AT level 
    (e.g. class). Two options: 
    - select an equal amount for each node at a given level.
    E.g. if we have 1000 class nodes, and we can afford to train for the total of
    10000 vectors, select 10 vectors for each class node. The problem with that is
    the potential overrepresentation of some branches of life above the class level 
    (those with more classes sequenced).
    How the sequence is collected for a given node: each node has a list of leaf
    nodes below it (LN). Suppose a node is allocated K vectors. For each i in 1...K
    pop a random element x from LN and select a random vector from x. If LN becomes
    empty, restore LN to its original value and continue till i==K. This procedure
    assures that we avoid sampling the same leaf twice while there are still unsampled 
    leaves, and at the same time that a random sample is picked up from available
    leaves. 
    - select an equal amount for each child at every node in AT. E.g. if a node is
     asigned 1000 vectors, and it has 5 children, then each child is asigned 200
     vectors. The same is repeated for each child - its allocated 200 vectors are
     split equally between its children. The problem with that: at some point we
     reach 1 vector allocated for the node that has several children. More general,
     we have N*l + k vectors allocated to the node, and N > k children, l can be
     zero. We select non-equal children for the k vectors at random. It is obvious
     that it does not make sense to train with a very few vectors for a given node
     when we have much more sequence available for it. So, we need to modify the
     selection algorithm: if the allocated number drops below certain cutoff level
     C for a given node in AT, we use C as the allocated number. Which brings 
     questions: why the tree structure above should influence the amount of sequence
     when we train for a certain level below? We probably need to assign equal amount
     to each class node when we train the class level classifier, and then go down
     in the manner described above in order to select vectors for each one.
     So, the modified procedure looks like this:
     Let us have N vectors total amount, for which we are able to train (e.g. 100000).
     Let K be the number of nodes at a given level X of AT (e.g. 1000 class level nodes)
     Set l = N/K to be the target allocated vectors for each node (1000 in our example).
     Starting with each X level node, recursively assign equal number of vectors to 
     each child of each node as described above. That procedure will construct the
     selection of vectors for each node at level X only. It will be repeated 
     independently for all other levels from AT. 

"""    


class NodeReducerVisitor:
    """Removes all nodes from a tree for which a given predicate is false.
    A node is removed by connecting all its children directly with its parent.
    It should be applied in the depth-top order.
    It never removes the starting node, so it must be already a valid reduced node
    (predicate is true) if you want to get the full subtree in a reduced form."""
    
    def __init__(self,predicate,tree):
        """@param predicate - if  not predicate(node), remove the node
        @param tree - the entire tree object; it will be used to clean up after node removal"""
        self.predicate = predicate
        self.tree = tree
    
    def __call__(self,node):
        # leaves will be removed while processing their parent
        if not node.isLeaf():
            predicate = self.predicate
            deleted = []
            #This will iterate by lifting grandchildren up until
            #all children of the current node are have 'predicate' to be true.
            #Then we can return and the recursion will proceed down to the 
            #newly formed children list.
            
            children_good = []
            children_bad = []
            
            for child in node.children:
                if predicate(child):
                    children_good.append(child)
                else:
                    children_bad.append(child)
            
            while len(children_bad) > 0:
                new_children_bad = []
                for child in children_bad:
                    deleted.append(child.id)
                    #assign children of child to this node itself and
                    #split them into good and bad lists
                    for grchild in child.children:
                        grchild.par = node
                        grchild.idpar = node.id
                        if predicate(grchild):
                            children_good.append(grchild)
                        else:
                            new_children_bad.append(grchild)
                children_bad = new_children_bad
            
            node.children = children_good
            self.tree.removeNodes(ids=deleted)


class RankFakerVisitor:
    """Assigns fake rank hierarchy starting from a top node, e.g. for uniformity of viral classification with other 'kingdoms'.
    Ranks are taken from a given list. The process goes from top to bottom and stops when leaf is found
    or if an already assignd 'species' rank is encountered. We want to preserve the originally assigned
    'species' ranks because the distinction between 'species' and underlying 'strains' can be important
    when we subsequently select genetic sequence for a given rank."""
    
    def __init__(self,topNode,lineage):
        self.topNode = topNode
        lineage = list(lineage)
        if lineage[-1] == 'species':
            lineage = lineage[:-1]
        self.lineage = lineage

    def __call__(self,node):
        #we go down to leaf nodes and reassign ranks of each leaf node lineage
        #this will repeatedly reassign ranks of internal nodes, but should be safe
        if node.isLeaf():
            lin_ranks_template = self.lineage
            lin = node.lineage(topNode=self.topNode)
            #make lineage top to bottom
            lin.reverse()
            lin_ranks = [ n.rank for n in lin ]
            #find position of 'species' rank in the lineage,
            #or set the last (this node) rank to be 'species'
            try:
                ind_species = lin_ranks.index('species')
            except ValueError:
                ind_species = len(lin_ranks) - 1
                lin_ranks[ind_species] = 'species'
            #set all other nodes in the lineage to 'no_rank'
            for i in range(ind_species):
                lin_ranks[i] = noRank
            for i in range(ind_species+1,len(lin_ranks)):
                lin_ranks[i] = noRank
            #set nodes up to the 'species' node to ranks from template lineage
            for i in range(min(ind_species,len(lin_ranks_template))):
                lin_ranks[i] = lin_ranks_template[i]
            #actually assign constructed ranks to nodes in the lineage
            for i in range(len(lin_ranks)):
                lin[i].rank = lin_ranks[i]
    
    def __call__old(self,node):
        lineage = list(self.lineage) # 'tuple' has no 'index' method
        par_rank = node.getParent().rank
        lin_orig = node.lineageRanksStr()
        rank = node.rank
        #do not do anything to the starting point,
        #just continue to children
        if rank == lineage[0]:
            pass
        #if the parent is 'species', then continue down and
        #set anything below to 'no_rank' (presumably strains)
        elif par_rank == 'species':
            if rank == 'species':
                print "Warning: 'species' is marked as a parent of 'species', original lineage  %s, current lineage %s" % \
                (lin_orig,node.lineageRanksStr())
            node.rank = noRank
        #do nothing if this node is already a 'species',
        #continue to the children
        elif rank == 'species':
            pass
        #only after we checked current for 'species', we can
        #propagate parents' 'no_rank'
        elif par_rank == noRank:
            node.rank = noRank
        else:
            next_rank_ind = lineage.index(par_rank) + 1
            if next_rank_ind < len(lineage):
                node.rank = lineage[next_rank_ind]
            else:
                node.rank = noRank
        #Set the rank to 'species' for every leaf node
        #that does not have a 'species' already in its lineage
        #ASSUMPTION: What if there are nodes in NCBI taxonomy tree
        #that are meant to be above species level but do not have
        #any children? We check here and warn if we are about to
        #reset explicit leaf rank (e.g. phylum), but we have no
        #clue about the meaning of 'no_rank' leaves. Perhaps other
        #fields in node.data can tell?
        if node.isLeaf():
            lin_before = node.lineageRanks()
            if 'species' not in lin_before:
                node.rank = 'species'
                if rank != noRank:
                    print ("Warning: resetting existing rank '%s' in leaf node to 'species' for taxid = %s,"+\
                           " original lineage: %s, current lineage: %s, previous lineage: %s") % \
                           (rank,node.id,lin_orig,node.lineageRanksStr(),lin_before)



class TaxaLevels:
    """Class that assigns classification levels to nodes of the taxonomic tree.
    We will train our classifiers to predict these levels."""

    unclassId = 10

    def __init__(self,taxaTree=None):
        """Constructor.
        @param taxaTree instance of TaxaTree - if not None, will be modified, so that TaxaLevels.lineage() etc work after words; 
        otherwise TaxaLevels.setLevels(taxaTree) should be called separately.
        """
        self.levels = list(linnMainRanks)
        self.levelIds = {noRank:0, unclassRank:self.unclassId, "species":20, "genus":30, "family":40,
                         "order":50, "class":60, "phylum":70, "kingdom":80, "superkingdom":90}
        self.levelByTaxid = { "superkingdom" : (viralRootTaxid,) }
        self.viralLevels = ("superkingdom","family","genus","species",unclassRank)
        self.levelSet = set(self.levels) | set((unclassRank,))
        #index of a given level name in the list 'levels'
        levelPos = {}
        for (i,level) in zip(range(len(self.levels)),self.levels):
            levelPos[level] = i
        self.levelPos = levelPos
        if taxaTree is not None:
            self.setLevels(taxaTree)

    def getLevelId(self,name):
        return self.levelIds[name]

    def getLevelNames(self,order="ascend"):
        if order == "ascend":
            return list(self.levels)
        elif order == "descend":
            l = list(self.levels)
            l.reverse()
            return l
        else:
            raise ValueError("Unknown value for 'order': " + str(order))

    def isNodeLevelLessOrEqual(self,node,level):
        return self.getLevelId(node.level) <= self.getLevelId(level)

    def setLevels(self,taxaTree):
        levelSet = self.levelSet
        for node in taxaTree.getNodesIter():
            if node.rank in levelSet:
                node.level = node.rank
            elif node.isUnclassified():
                node.level = unclassRank
            else:
                node.level = noRank
        for level in self.levelByTaxid:
            taxids = self.levelByTaxid[level]
            for taxid in taxids:
                node = taxaTree.getNode(taxid)
                node.level = level
        ##at this point, node.level is either one of 'levels' or 'no_rank'
        ##remove all viral levels above family (there are only two 'orders' currently defined
        ##for viruses, so it does not make much sense to train for them)
        viralLevels = self.viralLevels
        for node in taxaTree.iterDepthTop(viralRootTaxid):
            if node.level not in viralLevels:
                node.level = noRank
        levelIds = self.levelIds
        for node in taxaTree.getNodesIter():
            node.idlevel = levelIds[node.level]

    def lineage(self,node,withUnclass=True):
        levelSet = self.levelSet
        return [ n for n in node.lineage() if ( n.level in levelSet ) and (withUnclass or n.level != unclassRank)]

    def lineageKeys(self,node,getDistance=False):
        levelIds = self.levelIds
        nodeDepth = node.getDepth()
        if getDistance:
            return [ (levelIds[n.level],n.id,nodeDepth - n.getDepth()) for n in self.lineage(node)]
        else:
            return [ (levelIds[n.level],n.id) for n in self.lineage(node)]

    def lineageFixedList(self,node,null=None):
        """Return a list of taxids that correspond to the list of level names returned by getLevelNames("ascend").
        If a given level is not present in this node's lineage, the corresponding element is 'null'."""
        levelPos = self.levelPos
        ret = [null]*len(self.levels)
        for n in self.lineage(node,withUnclass=False):
            ret[levelPos[n.level]] = n.id
        return ret

    def reduceNodes(self,tree,topNodeId=None):
        tree.visitDepthTop(NodeReducerVisitor((lambda node: node.idlevel != 0),tree),topNodeId)
        tree.reindex()
        

def selectTaxaGis(taxids,taxaTree,giToTaxa,withSubnodes=True):
    tax = set(taxids)
    assert 0 not in tax
    if withSubnodes:
        for idnode in taxids:
            for node in taxaTree.iterDepthTop(idnode):
                tax.add(node.id)
    taxgi = []
    for gi in range(len(giToTaxa))[:10000]:
        taxid = giToTaxa[gi]
        if True or taxid in tax:
            taxgi.append((gi,taxid))
    return n.array(taxgi,dtype=[("gi","i8"),("taxid","i8")])

def downloadTaxaGenbankSeqs(taxids,outFile,taxaTree,giToTaxa,withSubnodes=True):
    from MGT.Entrez import EzRequest 
    taxgi = selectTaxaGis(taxids=taxids,
            taxaTree=taxaTree,
            giToTaxa=giToTaxa,
            withSubnodes=withSubnodes)
    pdb.set_trace()
    print "Getting Genbank records for %s GIs" % len(taxgi)
    gis = taxgi["gi"]
    req = EzRequest(batchSize=50000)
    outFile=open(outFile,'w')
    giChunk=40000
    for iGi in xrange(0,len(gis),giChunk):
        iGiEnd = min(iGi+giChunk,len(gis))
        req.fetch(ids=gis[iGi:iGiEnd],outFile=outFile)
        


