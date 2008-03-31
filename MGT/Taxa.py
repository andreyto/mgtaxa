"""Classes to represent and traverse a taxonomy tree (NCBI format)"""

from MGT.Common import *

import os
import numpy
from numpy import array
import itertools

class TaxaNodeData(object):
    fields = \
    (
        ('taxid',int),
        ('partaxid',int),
        ('rank',str),
#        ('embl_code',str),
#        ('divid',int),
#        ('inh_div',bool),
#        ('gcode_id',int),
#        ('inh_gc',bool),
#        ('mgcode_id',int),
#        ('inhmgc',bool),
#        ('gbhidden',bool),
#        ('hidsubtree',bool),
#        ('comments',str)
     )
    
    savedDtype = numpy.dtype([('taxid','int32'),
                              ('partaxid','int32'),
                              ('rank','S20')])
    
    def __init__(self,values):
        for (f,v) in zip(self.fields,values):
            setattr(self,f[0],f[1](v))
        self.rank = self.rank.replace(' ','_')
        
    def __str__(self):
        return ' | '.join([ str(getattr(self,f[0])) for f in self.fields ])
            
        
class TaxaNode(object):
    
    def __init__(self,values):
        self.data = TaxaNodeData(values)
        self.id = self.data.taxid
        self.children = []
        
    def __str__(self):
        return str(self.data)
        
    def setParent(self,par):
        #TODO: is garbage collector going to deal with circular refs?
        self.par = par
        if par is not None:
            par.children.append(self)
        
    def getParent(self):
        return self.par
    
    def getChildren(self):
        return self.children
        
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

    def lineageRanks(self,*l,**kw):
        return [ node.data.rank for node in self.lineage(*l,**kw) ]
    
    def lineageRanksTaxa(self,*l,**kw):
        return [ (node.data.rank,node.data.taxid) for node in self.lineage(*l,**kw) ]
    
    def lineageRanksStr(self,*l,**kw):
        return ','.join([ "%s=%s" % (rank,taxid) for (rank,taxid) in self.lineageRanksTaxa(*l,**kw) ])

    def visitDepthTop(self,func):
        """If functor 'func' does not return anything when applied to itself, recursively do the same with each of the children."""
        if func(self) is None:
            for child in self.children:
                child.visitDepthTop(func)
    
    def iterDepthTop(self):
        """Return depth-first top-to-bottom iterator."""
        yield self
        for child in self.children:
            for node in child.iterDepthTop():
                yield node


    def visitDepthBottom(self,func):
        """If functor 'func' does not return anything when applied to itself, recursively do the same with each of the children."""
        for child in self.children:
            child.visitDepthBottom(func)
        func(self)


class TaxaTree(object):
    
    def __init__(self,nodes=None,ncbiDumpFile=None,save=False,load=False):
        #It is still much faster to load node file directly rather than
        #from a pickled version of raw (unlinked) nodes or from a pickled
        #numpy array, and uses much less
        #memory. Therefore, the default values for cacheing are False.
        assert nodes is None or ncbiDumpFile is None
        assert not (ncbiDumpFile is None and (save or load))
        #pickling the nodes after they were internally linked in a tree structure
        #leads to huge memory consumption and very slow.
        
        if ncbiDumpFile is not None:
            cacheFile = os.path.basename(ncbiDumpFile)+'.cache'
            if load and os.path.isfile(cacheFile):
                print "Loading cached version of TaxaTree"
                nodes = self.loadData(cacheFile)
            else:
                inp = open(ncbiDumpFile,'r')
                nodes = {}
                n_splits = len(TaxaNodeData.fields)
                for rec in inp:
                    node = TaxaNode([ x.strip() for x in rec.split('\t|\t',n_splits)[:n_splits]])
                    nodes[node.id] = node
                inp.close()
                if save:
                    self.saveData(nodes,cacheFile)
            
        self.nodes = nodes
        self.rootNode = None
        for v in nodes.itervalues():
            partaxid = v.data.partaxid
            #in NCBI file, root node points to itself as a parent
            #we use None
            if partaxid != v.data.taxid:
                par = nodes[partaxid]
            else:
                par = None
                self.rootNode = v
            v.setParent(par)
    
    def saveData(self,nodes,fileName):
        arr = numpy.empty(len(nodes),dtype=TaxaNodeData.savedDtype)
        i = 0
        for node in nodes.itervalues():
            arr[i] = (node.data.taxid,node.data.partaxid,node.data.rank)
            i += 1
        dumpObj(arr,fileName)
    
    def loadData(self,fileName):
        arr = loadObj(fileName)
        nodes = {}
        for rec in arr:
            nodes[int(rec["taxid"])] = TaxaNode(rec)
        return nodes

    def write(self,out):
        for id in sorted(self.nodes.iterkeys()):
             out.write(str(self.nodes[id])+'\n')
    
    def writeLineage(self,out):
        for id in sorted(self.nodes.iterkeys()):        
            out.write(self.nodes[id].lineageRanksStr()+'\n')
             
    def getNode(self,id):
        return self.nodes[id]
    
    def getNodesDict(self):
        return self.nodes
    
    def getNodesIter(self):
        return self.nodes.itervalues()
    
    def getRootNode(self):
        return self.rootNode

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

        
    def visitDepthBottom(self,func,id=None):
        if id is None:
            top = self.getRootNode()
        else:
            top = self.getNode(id)
        top.visitDepthBottom(func)

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
            node_rank = node.data.rank
            if node_rank == rank:
                mask[node_id] = node_id
            elif node.par is not None:
                mask[node_id] = mask[node.par.id]
        self.visitDepthTop(mask_setter)
        return mask


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

class RankReducerVisitor:
    """Removes all internal (non-leaf, non-root) nodes from a tree with rank names not in a given list.
    Example of such list is a 'linnMainRanks' module level variable ('domain',...,'species').
    A node is removed by connecting all its children directly with its parent.
    It should be applied in the depth-top order.
    It does not check the starting node, so it must be already valid reduced node if you
    want to get the full subtree in a reduced form. Root node is a valid starting point."""
    
    def __init__(self,allowedRanks):
        self.allowedRanks = allowedRanks
    
    def __call__(self,node):
        if not node.isLeaf():
            allowedRanks = self.allowedRanks
            
            #This will iterate by lifting grandchildren up until
            #all children of the current node are allowedRanks or leaves.
            #Then we can return and the recursion will proceed down to the 
            #newly formed children list.
            
            children_good = []
            children_bad = []
            
            for child in node.children:
                if child.isLeaf() or child.data.rank in allowedRanks:
                    children_good.append(child)
                else:
                    children_bad.append(child)
            
            while len(children_bad) > 0:
                new_children_bad = []
                for child in children_bad:
                    #assign children of child to this node itself and
                    #split them into good and bad lists
                    for grchild in child.children:
                        grchild.par = node
                        if grchild.isLeaf() or grchild.data.rank in allowedRanks:
                            children_good.append(grchild)
                        else:
                            new_children_bad.append(grchild)
                children_bad = new_children_bad
            
            node.children = children_good

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
            lin_ranks = [ n.data.rank for n in lin ]
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
                lin[i].data.rank = lin_ranks[i]
    
    def __call__old(self,node):
        lineage = list(self.lineage) # 'tuple' has no 'index' method
        par_rank = node.getParent().data.rank
        lin_orig = node.lineageRanksStr()
        rank = node.data.rank
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
            node.data.rank = noRank
        #do nothing if this node is already a 'species',
        #continue to the children
        elif rank == 'species':
            pass
        #only after we checked current for 'species', we can
        #propagate parents' 'no_rank'
        elif par_rank == noRank:
            node.data.rank = noRank
        else:
            next_rank_ind = lineage.index(par_rank) + 1
            if next_rank_ind < len(lineage):
                node.data.rank = lineage[next_rank_ind]
            else:
                node.data.rank = noRank
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
                node.data.rank = 'species'
                if rank != noRank:
                    print ("Warning: resetting existing rank '%s' in leaf node to 'species' for taxid = %s,"+\
                           " original lineage: %s, current lineage: %s, previous lineage: %s") % \
                           (rank,node.id,lin_orig,node.lineageRanksStr(),lin_before)


viralRootTaxid = 10239
#main Linnaean ranks, modified with superkingdom in place of domain rank
linnMainRanks = ("species","genus","family","order","class","phylum","superkingdom")
viralRanksTemplate = linnMainRanks
noRank = "no_rank"

class TaxaLevels:
    """Class that assigns classification levels to nodes of the taxonomic tree.
    We will train our classifiers to predict these levels."""

    def __init__(self):
        self.levels = list(linnMainRanks)
        self.levelIds = {noRank:0, "species":20, "genus":30, "family":40,
                         "order":50, "class":60, "phylum":70, "superkingdom":80}
        self.levelByTaxid = { "superkingdom" : (viralRootTaxid,) }
        self.viralLevels = ("superkingdom","family","genus","species")
        self.levelSet = set(self.levels)
        #index of a given level name in the list 'levels'
        levelPos = {}
        for (i,level) in zip(range(len(self.levels)),self.levels):
            levelPos[level] = i
        self.levelPos = levelPos

    def getLevelNames(self):
        return self.levels

    def setLevels(self,taxaTree):
        levelSet = self.levelSet
        for node in taxaTree.getNodesIter():
            if node.data.rank in levelSet:
                node.data.level = node.data.rank
            else:
                node.data.level = noRank
        for level in self.levelByTaxid:
            taxids = self.levelByTaxid[level]
            for taxid in taxids:
                node = taxaTree.getNode(taxid)
                node.data.level = level
        ##at this point, node.data.level is either one of 'levels' or 'no_rank'
        ##remove all viral levels above family (there are only two 'orders' currently defined
        ##for viruses, so it does not make much sense to train for them)
        viralLevels = self.viralLevels
        for node in taxaTree.iterDepthTop(viralRootTaxid):
            if node.data.level not in viralLevels:
                node.data.level = noRank

    def lineage(self,node):
        levelSet = self.levelSet
        return [ n for n in node.lineage() if n.data.level in levelSet ]

    def lineageKeys(self,node):
        levelIds = self.levelIds
        return [ (levelIds[n.data.level],n.id) for n in self.lineage(node)]

    def lineageFixedList(self,node):
        """Return a list of taxids that correspond to the list of level names returned by getLevelNames().
        If a given level is not present in this node's lineage, the corresponding element is None."""
        levelPos = self.levelPos
        ret = [None]*len(self.levels)
        for n in self.lineage(node):
            ret[levelPos[n.data.level]] = n.id
        return ret