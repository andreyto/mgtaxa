"""Classes that assign classification labels to various ranks in the taxonomy of our known sequence.
We adopt a hierarchical approach to classification in order to place a sequence onto a known
taxonomic tree. A recursive procedure is applied to each sample: if a node at some rank has
been assigned to it at a previous step, then a classifier trained for only sub-nodes of this
node at the next lower rank is applied to assign a sub-node. E.g. if a node is assigned
a given order, we apply a classifier trained to assign families under this order.
In addition, a one-class classifier is trained for each sub-node (family in our example).
If its prediction disagrees with the multi-class prediction, the sample is assigned to the
"reject" group.
The code in this module:
    - groups available training sequence under each node of each classification ranks
    - assigns classification labels based on ranks with addition of special label for
    sequences that do not have a required rank assigned to them (e.g. a sequence has
    order and genus but no family aasigned)
    - splits sequence into balanced testing and training sets, restricted by maximum
    and minimum sequence lengths and excluding from training all genera selected for testing
    - the testing sequence length in turn reflects the expected length distribution
    of prediction sequence fragments
    - after the testing sequence is set aside, the splitting procedure is again applied
    to training sequence in order to create datasets for n-fold parameter optimization.

Final output:
A set of tables that describe a hierarchy of classifiers and their training/testing sequence.
Structure of these tables allow straightforward generation of classification training/validation
jobs for each classifier independently, as well as testing/prediction with the entire hierarchical
algorithm.

Implementation:
We describe first a set of operations and then determine how to map them into
SQL or in-memory tree representation.
Terminology: in order to avoid confusion between taxonomic 'class' rank and
classification 'class', we replace the later with 'label' in the following text.

1. For each node at a given level (rank from our predefined set of classification
taxonomic ranks), select all next lower level sub-nodes (e.g. all 'family'
nodes under a given 'order' node). Direct subnodes with missing
next lower level rank form one separate group "unknown". Each group is assigned its
taxid as a class label, except for "unknown" group which is assigned label 1. 1 is an actual
taxid of NCBI root node, but that node will never appear as a classification target,
so it is safe to re-use it for arbitrary label. Label for the "other" node should be
made a parameter just in case.
Special case: it is possible that there are no subnodes with the next lower level
rank (e.g. we eliminated the only two existing 'order' nodes from the viral
hierarchy in order to deal with 'family' directly under viral 'superkingdom'.
And there are no 'phylum' or 'class' nodes for viruses.). The implementation
should use the next-next-lower level rank and proceed.
Special case: only one sub-node. We should replace it with its next level sub-nodes.

Discussion on the meaning of "unknown" group:
Two diferent origins are possible:
    a) The "real" classification hierarchy does not actually have a given rank defined
    for these sequences. Example are viral 'order' (only two orders are defined). Although
    we dealt with the viral case at the earlier stage by eliminating the "order" rank
    completely, such a case can be present at other branches of the tree. In that case
    it appears to be OK for us to lump all such sequences in a separate group, train
    to classify it along with "real" labels, and then apply the standard recursive procedure
    to that group. For example, if all viral families are separable from each other, then
    any non-intersecting sets of them must be also separable, and thus it is safe to
    classify first at the 'order' level (two defined orders and one "unknown" order
    with the rest of the families), and then at the 'family' level under each of the
    three 'order' labels.
    b) The required rank was not assigned to a given subnode for some reason (e.g. lack
    of information about its phenotype), but in reality it belongs to some of the defined
    nodes at this rank. E.g. for a given family node, no 'order' was defined (it is attached
    directly to some 'class' node), but given full information about the constituent
    species, it would be attached to one of the existing 'order' nodes. In such case creating
    a separate training label is wrong - such label will overlap with other labels and lead
    to misclassifications (given our central assumption that sequence composition reflects
    taxonomic hierarchy).
"""


from MGT.Common import *
from MGT.Taxa import *
from MGT.Sql import *
from MGT.SampDb import *
from MGT.SampDbKmer import *
from MGT.PredictorDb import *

def exclude_rRNA(self):
    db.createTableAs(
    """select * from seq_hdr where hdr rlike '[[:<:]]rRNA[[:>:]]'"""
    )


class SamplerConcat(MGTOptions):
    """Sample chunks of available training sequence by concatenating all sequences for each taxid.
    Concatenating (as opposed to chunking each sequence individually) greatly simplifies random
    sampling of taxid chunks w/o replacement and increases the amount of sample data for small
    sequences. It creates chimeric k-mer vectors however (but not chimeric k-mers because we insert
    spacers)."""

    # Subtrees starting at these nodes are considered "testing terminal" -
    # entire subtree is always selected for testing
    # We need "sub..." ranks because there subspecies that do not have
    # a species or genus node in their lineage.

    termTopRanks = ("genus","species","subspecies","subgenus")

    def __init__(self,db,sampLen,kmerLen,namePrefix):
        """Constructor.
        @param db - SqlDB instance
        @param sampLen - length of sequence chunks to sample
        @param kmerLen - length of k-mer (k)
        @param namePrefix - files will be in a directory with this name and tables created by this 
        object will have this prefix in their names.
        """
        MGTOptions.__init__(self)
        self.stage = 2
        self.db = db
        self.sampLen = sampLen
        self.kmerLen = kmerLen
        self.namePrefix = namePrefix
        self.tblPrefix = namePrefix + '_'
        self.dirPrefix = namePrefix
        makedir(self.dirPrefix)
        self.hdfSampFile = self.getFileName("samp.hdf")
        # @todo: after we debug this faster version, load the tree from SQL storage to reflect possible
        # tree modifications saved from earlier stages
        self.storePickle = NodeStoragePickleSep(fileName=self.getFileName(name="taxaTree.postMarkTraining.pkl"))
        self.storeDump = NodeStorageNcbiDump(ncbiDumpFile=self.taxaNodesFile,
                ncbiNamesDumpFile=self.taxaNamesFile)
        self.levels = TaxaLevels()
        if self.stage == 1:
            print "DEBUG: Loading the tree from dump file"
            self.taxaTree = TaxaTree(storage=self.storeDump)
            self.levels.setLevels(self.taxaTree)
        else:
            print "DEBUG: Loading the tree from pickle"
            self.taxaTree = TaxaTree(storage=self.storePickle)
        print "DEBUG: tree loaded"
        self.taxaTreeStore = NodeStorageDb(db=self.db,tableSfx=self.taxaTreeTableSfxMain)

    def getTableName(self,name):
        """Given a stem table name, return full table name unique for this object (prefix+name).
        @param name - stem name
        """
        return self.tblPrefix+name.strip()

    def getFileName(self,name):
        """Same as getTableName() but for files."""
        return os.path.join(self.dirPrefix,name)

    def rebuild(self):
        if self.stage == 1:
            self.mkHdfSampleInd()
            self.mkDbSampleCounts()
            print "DEBUG: loadSampleCountsMem()"
            self.loadSampleCountsMem()
            print "DEBUG: setIsUnderUnclassMem()"
            self.setIsUnderUnclassMem()
            print "DEBUG: markTrainable()"
            self.markTrainable()
            print "DEBUG: selectTestTaxa()"
            self.selectTestTaxa()
            print "DEBUG: markTraining()"
            self.markTraining()
            print "DEBUG: checkPostMarkTraining()"
            self.checkPostMarkTraining()
            print "DEBUG: saving the tree"
            self.taxaTree.delAttribute('sampSel')
            #raise ValueError(0)
            self.storePickle.save(tree=self.taxaTree)
            print "DEBUG: Writing statistics into SQL"
            self.statTestTaxa()
        else:
            print "DEBUG: writeTraining()"
            #self.tmp_hdfSamplesConv()
            #return
            self.setTrainTarget()
            trNodes = self.selectMicFamilies()
            self.writeTraining(trNodes=trNodes)
            self.writeTesting()
            self.tmp_hdfSamplesConv()

    def mkHdfSampleInd(self):
        hdfFileActSeq = pt.openFile(self.hdfActSeqFile,mode="r")
        hdfSeqInd = hdfFileActSeq.getNode(self.hdfActSeqInd)
        hdfFileSamp = pt.openFile(self.hdfSampFile,mode="w")
        hdfMakeSampIndexConcat(hdfFileSamp,self.hdfSampGroup,hdfSeqInd,self.sampLen)
        hdfFileSamp.close()
        hdfFileActSeq.close()

    def mkDbSampleCounts(self):
        """Make DB tables with sample counts per taxonomy id"""
        db = self.db
        # Because we concatenate sequence, we can get sample counts
        # from act_src
        tblTxSamp = self.getTableName("tx_samp")
        db.createTableAs(tblTxSamp,"""
        select taxid,floor(sum(seq_len)/%i) as samp_n from act_src group by taxid having samp_n > 0
        """ % self.sampLen)
        db.createIndices(primary="taxid",table=tblTxSamp)
        db.analyze(tblTxSamp)

    def loadSampleCountsMem(self):
        """Transfer sample counts from DB into in-memory taxonomy tree object."""
        tblTxSamp = self.getTableName("tx_samp")
        loadDb = True
        self.taxaTree.getRootNode().setIndex()
        if loadDb:
            self.taxaTreeStore.loadAttribute(tree=self.taxaTree,
                                         name="samp_n",
                                         sql="select taxid,samp_n from %s" % tblTxSamp,
                                         default=0L,
                                         setDefault=True,
                                         ignoreKeyError=True,
                                         typeCast=long)

            samp_n = self.taxaTree.getRootNode().makePropArrayFromAttrib(dtype='i4',name='samp_n')
            dumpObj(samp_n,"samp_n.pkl")
        else:
            samp_n = loadObj("samp_n.pkl")
            print "DEBUG: samp_n.shape = ", samp_n.shape
            self.taxaTree.getRootNode().setAttribFromPropArray(prop=samp_n,name="samp_n")
        self.taxaTree.setTotal("samp_n","samp_n_tot")
        # samp_n_tot_class excludes immediate unclassified children but includes
        # unclassified grand-children
        for node in self.taxaTree.iterDepthTop():
            node.samp_n_tot_class = node.samp_n_tot - \
                    sum( ( child.samp_n_tot for child in node.getChildren() if child.isUnclassified() ), 0 )
        #self.taxaTree.setReduction("samp_n","samp_n_tot_class",condition=lambda node: not node.isUnclassified()) 
    
    def setIsUnderUnclassMem(self):
        """Set an attribute 'isUnderUnclass' for in-memory tree nodes.
        It flags all nodes that have "unclassified" super-node somewhere in their lineage."""
        taxaTree = self.taxaTree
        iter = taxaTree.iterDepthTop()
        root = iter.next()
        root.isUnderUnclass = False
        for node in iter:
            if node.isUnclassified() or node.getParent().isUnderUnclass:
                node.isUnderUnclass = True
            else:
                node.isUnderUnclass = False
        
    def markTrainable(self):
        """Mark all nodes that have enough samples for training with and without validation.
        When we train for subsequent validation with testing set, we cannot use any 'unclassified'
        samples because they might actually represent the excluded testing set. On the other
        hand, there is no reason not to use any 'unclassified' samples under a 'classified'
        node for the final re-training.
        @pre loadSampleCountsMem() and setInUnderUnclassMem() were called
        @post in-memory taxaTree nodes have the following attributes set:
            isTrainable - trainable with 'unclassifed' samples included
            isTrainableTest - trainable for validation. A node can still be excluded from
            testing later if extracting testing set will leave too few samples for training,
            but it will be trained during testing stage - just will not have any true positive 
            testing samples."""
       
        self.setSampSel()
        taxaTree = self.taxaTree
        for node in taxaTree.iterDepthTop():
            # a quick filter that excludes nodes w/o any samples
            # nodes under 'unclassfied' subtrees are never the training targets
            if node.samp_n_tot > 0 and not node.isUnderUnclass:
                node.isTrainable = ( node.samp_n_tot >= node.sampSel.prod.train.min )
                testSampSel = node.sampSel.test
                node.isTrainableTest = ( node.samp_n_tot_class >= (testSampSel.train.min + testSampSel.test.min) )
                testTarget = max(\
                        min(\
                            int(round(node.samp_n_tot_class * testSampSel.test.ratio)),
                            node.samp_n_tot_class - testSampSel.train.min\
                            ),
                        testSampSel.test.min)
                node.nSampTestMax = testTarget
            else:
                node.isTrainable = False
                node.isTrainableTest = False
     
    def setSampSel(self):
        taxaTree = self.taxaTree
        viralRoot = taxaTree.getNode(viralRootTaxid)
        sampSel = self.sampSel
        for node in taxaTree.iterDepthTop():
            
            if node.isSubnode(viralRoot):
                node.sampSel = sampSel.vir
            else:
                node.sampSel = sampSel.all

    def mkDbTestingTaxa(self,rank="genus",superRank="family"):
        """Select complete 'rank' nodes uniformly distributed across 'superRank' super-nodes."""
        pairs = self.lowRankTestPairs(rank=rank,superRank=superRank)
        self.sampleTestPairs(pairs)

    def lowRankTestPairs(self,rank,superRank):
        """Return set of unique pairs (testing node,first supernode).
        @param rank - entire nodes of this rank will be selected as testing samples (e.g. genus).
        @param superRank - a second member of each pair will have at least this rank (e.g. family).
        'superRank' nodes might not be 'isTrainable'. That should be filtered later."""

        rankId = self.levels.getLevelId(rank)
        superRankId = self.levels.getLevelId(superRank)
        noRankId = self.levels.getLevelId(noRank)
        #select all nodes with samples and their superrank super-nodes
        #@todo idlevel is more coarse than it has to be for this purpose.
        #ideally, we should pick superfamily if we could not find family -
        # with current idlevels we will pick order.
        taxaTree = self.taxaTree
        pairs = set()
        for node in taxaTree.iterDepthTop():
            if node.samp_n > 0 and not node.isUnderUnclass:
                sub = node
                sup = None
                for n in node.lineage():
                    if n.idlevel != noRankId:
                        #nodes as high as rankId will be selected in entirety
                        if n.idlevel <= rankId:
                            sub = n
                        #nodes at least as high as superRankId will be uniformly sampled
                        if n.idlevel >= superRankId:
                            sup = n
                            break
                # this should always be true, but just in case
                if sub.isSubnode(sup) and sup is not None:
                    #adding pairs into a set object makes sure
                    #that they are unique (we can hit one genus many times ascending from its species)
                    pairs.add((sub,sup))
        db = self.db
        db.saveRecords(records = ( (sub.id,sup.id) for (sub,sup) in pairs ), 
                       table = SqlTable(name="test_pairs",
                                        fields = (SqlField("id_sub","integer"),SqlField("id_sup","integer"))))
        return pairs

    def sampleTestPairs(self,pairs):
        """Randomly sample a list of (testing node,first supernode) pairs obtained with lowRankTestPairs()."""
        #samp_n_tot - test check above for every genus
        pass

    def selectTestTaxa(self):
        """Select nodes that will be used entirely for testing.
        The method aims to achieve these goals:
        - create a uniform random sampling of the tree
        - sample as deep into the hierarchy as possible
        - select entire genera or entire species if there is no genera above them
        (we call such nodes "terminal" here)
        This is done recursively: for a given node: 
        a) all taxa selected as testing for its non-terminal children is also selected
        b) and then terminal children are added at random until we reach the limit on 
        testing samples set for this node.
        The major alternative would be to do the (b) selection from a list of all
        so far non-selected terminal subnodes in an entire sub-tree of a given node
        (the method would have to return such a list from each recursive call, with the
        parent concatenating the lists). This would provided better chances to reach the
        testing ratio target for higher order taxa. However, this would also run a danger
        of a non-uniform training and deteriorating perfomance. For example, if we have
        a subtree (A,(B,C),D) and B was selected as testing for (B,C) subtree, we could also
        select C as testing for (A,(B,C),D) subtree (it would still be used as training
        for (B,C). However, we would have no trainig samples from (B,C) when trainig (A,(B,C),D),
        which can become an issue.
        @post Attributes are set for each node:
            - is_test - if True, entire subtree is selected as testing sample. Set for every node
            in such a subtree.
            - has_test - if True, node has enough testing samples, but it is not 'is_test' itself.
            - samp_n_tot_test - number of testing samples for this node
            - samp_n_tot_train - number of training samples for this node.
        Note that the presence of "unclassified" sequence still creates problems:
        selection of even one testing sample excludes from training all unclassified nodes
        which are immediate children of all nodes in testing sample lineage. Thus, any of these nodes
        can become non-trainable even if it is trainable w/o testing. However, we need to accumulate
        the testing sequence for use by upper nodes, so this seems unavoidable. To check the effect
        of this tradeoff, temporarily make TaxaNode.isUnclassified() to always return False and re-run
        the testing/training selection. So far this experiment shown little impact."""

        taxaTree = self.taxaTree
        rootNode = taxaTree.getRootNode()
        for node in taxaTree.iterDepthTop():
            node.is_test = False
            node.has_test = False
            node.samp_n_tot_test = 0
            if not node.isUnderUnclass:
                node.samp_n_tot_train = node.samp_n_tot
            else:
                node.samp_n_tot_train = 0

        def actor(node):
            # if a node is not even trainable in testing phase with minimum selection of testing samples, 
            # set the number of testing samples to zero and return
            if node.isUnderUnclass or node.samp_n_tot == 0:
            #if node.samp_n_tot == 0:
                node.samp_n_tot_test = 0
                return
            # Collect all classified child nodes with non-zero classified sample counts, for convenience
            # This will be retained on stack during recursion, but our tree is not deep
            children = [ child for child in node.getChildren() if not child.isUnclassified() and child.samp_n_tot > 0 ] 
            # Genus and species subnodes, we mark as "testing terminal" meaning that we
            # can only select the entire subnode for testing (we could actually only exclude
            # from training entire subnodes, but allocate only part of their samples for testing,
            # which would result in spreading the testing sample over wider set of subnodes. We
            # do not do this, however, because that would reduce the feature space covered by
            # training samples).
            # For other subnodes, we call the same function recursively
            for child in children:
                if child.rank in self.termTopRanks:
                    child.isTestTerm = True
                    child.samp_n_tot_test = 0
                    child.samp_n_tot_train = child.samp_n_tot
                else:
                    assert child.samp_n == 0, \
                            "We currently cannot tolerate nodes that have directly attached sequence "+\
                            "and are not 'testing terminal'. Examples would be some sequence assigned "+\
                            "a family taxid rather than taxid of some species or 'unclassified' "+\
                            "subnode under that family."
                    child.isTestTerm = False
                    actor(child)
                    # there will be cases when it is not possible to select any testing samples 
                    # from a subnode within provided constrains, and the total number of samples
                    # is more than zero but less than the minimal training number.
                    # In such case, we mark it as "testing terminal", free for selection,
                    # 
                    if child.samp_n_tot_test == 0 and child.samp_n_tot_train < child.sampSel.test.train.min:
                        child.isTestTerm = True
            termChildren = []
            samp_n_tot_test = 0
            samp_n_tot_train = 0
            for child in children:
                if child.isTestTerm:
                    termChildren.append(child)
                samp_n_tot_test += child.samp_n_tot_test
                # Initially we set node's training to be a sum of traning of all subnodes,
                # including all training terminal ones. Then we will be subtracting the counts
                # of training terminals when we select them for testing
                samp_n_tot_train += child.samp_n_tot_train
            #assert node.samp_n_tot_class >= samp_n_tot_test,
            #assert node.samp_n_tot_class - samp_n_tot_test >= node.sampSel.test.train.min, \
            #    "Total number of samples selected for subnodes should never exceed"+\
            #    " the total number of classified samples for the parent minus min training size."
            # we selected samp_n_tot_test samples from subnodes so far
            # the node will be considered testable if it has at least node.nSampTestMin test samples,
            # with the target number of node.nSampTestMax and node.nNodesTestMax (which are computed earlier as a fixed
            # percentage of all classified samples under the node but bound by the min number of training nodes to be left).
            # Here, we sample randomly w/o substitution the "testing terminal" subnodes until we reach either one of the 
            # upper bounds or exhaust the subnodes list.
            nSampTestToAdd = node.nSampTestMax - samp_n_tot_test
            # That assert is no longer applicable after we added code around 'new_test' below - now we allow
            # to accumulate the amount of testing samples above the target ratio if this is the only way to
            # select at least the minimum amount of testing (assuming that the remaining training is still
            # sufficient)
            # assert nSampTestToAdd >= 0, "Testing sample selection inside subnodes should not exceed parent's limit"
            if nSampTestToAdd > 0:
                nrnd.shuffle(termChildren)
                #termChildren.sort(key=lambda n: n.samp_n_tot)
                for child in termChildren:
                    # We continue through the remaining children if adding the current one 
                    # would exceed the selection limit. That means, the procedure is biased
                    # somewhat toward the smaller nodes and not uniformly random anymore. This seems
                    # to be a fair traidoff - otherwise we would have smaller number of testable nodes.
                    # It is better to have some testing than no testing.
                    if child.samp_n_tot <= nSampTestToAdd and \
                        samp_n_tot_train - child.samp_n_tot >= node.sampSel.test.train.min:
                        # When adding whole nodes for testing, we should not ignore "unclassified" content under them
                        nSampTestToAdd -= child.samp_n_tot
                        samp_n_tot_test += child.samp_n_tot
                        samp_n_tot_train -= child.samp_n_tot
                        # We currently just mark the child as selected for testing, and then
                        # collect all such in another pass through the tree.
                        # Calling the logging method here would be a more flexible solution because
                        # it would not assume that the testing node is a direct child of a tested node, e.g.
                        # self.logTestSamples(node,child)
                        # But it would require more extensive coding for itself and other parts of the 
                        # algorithm such as selection of the training set.
                        for n in child.iterDepthTop():
                            n.is_test = True
                            n.has_test = False
                            n.samp_n_tot_test = n.samp_n_tot
                            n.samp_n_tot_train = 0
                # A node might have two (or more) terminal subnodes, each larger than
                # the target testing size. In that case we just take the smaller one
                # assuming that remaining training size is above minimum
                if samp_n_tot_test < node.sampSel.test.test.min:
                    remTermChildren = [ child for child in termChildren if not child.is_test ]
                    if len(remTermChildren) > 0:
                        remTermChildren.sort(key=lambda n: n.samp_n_tot)
                        child = remTermChildren[0]
                        new_test = samp_n_tot_test + child.samp_n_tot
                        new_train = samp_n_tot_train - child.samp_n_tot 
                        if new_test >= node.sampSel.test.test.min and \
                            new_train >= node.sampSel.test.train.min:
                            samp_n_tot_test = new_test
                            samp_n_tot_train = new_train
                            for n in child.iterDepthTop():
                                n.is_test = True
                                n.has_test = False
                                n.samp_n_tot_test = n.samp_n_tot
                                n.samp_n_tot_train = 0



            node.samp_n_tot_test = samp_n_tot_test
            node.has_test = ( samp_n_tot_test >= node.sampSel.test.test.min )
            node.samp_n_tot_train = samp_n_tot_train
            # If no testing sequence was selected, we can use all (with "unclassified")
            # sequence for training
            if samp_n_tot_test == 0 and not node.isUnderUnclass:
                for n in node.iterDepthTop():
                    n.samp_n_tot_train = n.samp_n_tot


        def debugOnNodeUpdate(node,name,value):
            if node.id == 10781 or node.idpar == 10781:
                print node.id, name, value

        #taxaTree.setDebugOnUpdate(debugOnNodeUpdate)
        actor(rootNode)
        iter = taxaTree.iterDepthTop()
        iter.next() # skip the root
        # if there is no testing sequence and it's not unclassified or it is unclassified
        # under a classified node with no testing sequence, use all samples for training,
        # including unclassified
        for node in iter:
            if node.samp_n_tot_test == 0 and \
                    (not node.isUnderUnclass or \
                    (node.getParent().samp_n_tot_test == 0 and \
                    node.getParent().samp_n_tot_train != 0)):
                node.samp_n_tot_train = node.samp_n_tot
        #taxaTree.unsetDebugOnUpdate()


    def markTraining(self):
        """Mark final training state for each node - 'is_class' and 'n_mclass' attributes.
        Nodes that are not 'is_class' are ignored during training.
        'is_class' means that a node can be used as a label in training the
        multiclass classifier for the parent node.
        'n_mclass' is also set to the number of classes below the node.
        'samp_n_max' is the maximum leaf sample count under this node.
        This method must be called after selectTestTaxa() and uses node attributes set by
        that method."""

        taxaTree = self.taxaTree
        iter = taxaTree.iterDepthTop()
        rootNode = iter.next()
        rootNode.is_class = True # need this for the condition below to work
        for node in iter:
            node.is_class = False
            par = node.getParent()
            # that will exclude every node below genus or species
            if par.is_class and not par.rank in self.termTopRanks:
                if node.samp_n_tot_train >= node.sampSel.test.train.min and \
                        not node.isUnderUnclass:
                    node.is_class = True
        for node in taxaTree.iterDepthTop():
            node.n_mclass = len([child for child in node.getChildren() if child.is_class])
        taxaTree.setReduction("samp_n",
                "samp_n_max",
                reduction=max,
                condition=lambda node: not node.isUnclassified()) 


    def setTrainTarget(self):
        self.setSampSel()
        taxaTree = self.taxaTree 
        sampLen = self.sampLen
        for node in taxaTree.iterDepthTop():
            if node.n_mclass > 1:
                child_targ = int(node.sampSel.test.train.max / node.n_mclass)
                for child in node.getChildren():
                    child.samp_n_train_targ = min(child_targ,child.samp_n_tot_train)
            node.has_long_seq = node.samp_n_max * sampLen >= node.sampSel.longSeq
            node.has_long_seq_self = node.samp_n * sampLen >= node.sampSel.longSeq

    def checkPostMarkTraining(self):
        """Consistency check after a call to markTraining()."""
        taxaTree = self.taxaTree
        iter = taxaTree.iterDepthTop()
        rootNode = iter.next()
        for node in iter:
            par = node.getParent()
            assert par.samp_n_tot_train >= node.samp_n_tot_train
            assert not (not par.is_class and node.is_class)
            assert not (node.samp_n > 0 and node.has_test), \
                    "Node with directly attached sequence cannot be marked as testable"
            assert node.samp_n != 0 or \
                    node.samp_n_tot_train == sum([ n.samp_n_tot_train for n in node.getChildren() ])

    def statTestTaxa(self):
        """Save SQL statistics about testing/training set selection.
        Call after markTraining()"""
        tblStat = self.getTableName("st_test")
        taxaTree = self.taxaTree
        self.taxaTreeStore.saveAttributes(tree=taxaTree,
            nameToSql={'samp_n_tot_test':{'type':'integer'},
                       'samp_n_tot_class':{'type':'integer'},
                       'samp_n_tot_train':{'type':'integer'},
                       'samp_n':{'type':'integer'},
                       'samp_n_max':{'type':'integer'},
                       'isTrainableTest':{'type':'bool','name':'is_trainable_test'},
                       'is_test':{'type':'bool'},
                       'has_test':{'type':'bool'},
                       'is_class':{'type':'bool'},
                       'n_mclass':{'type':'integer'},
                       'isTestTerm':{'type':'bool','name':'is_test_term'},
                       'nSampTestMax':{'type':'integer','name':'n_samp_test_max'},
                       'isUnderUnclass':{'type':'bool','name':'is_under_unclass'},
                       'lind':{'type':'integer'},
                       'rind':{'type':'integer'},
                       'lnest':{'type':'integer'},
                       'rnest':{'type':'integer'},
                       'depth':{'type':'integer'},
                       'rank':{'type':'char(20)'},
                       'idpar':{'type':'integer'}},
            table=tblStat)
        storeTreeView = NodeStorageHypView(fileName=tblStat+".hv3",
            labeler=lambda n: "%s_%s_%s_N_%s_R_%s_T_%s_T_%s" % (n.id,n.rank,n.name[:10],n.numChildren(),
                n.samp_n_tot_train,n.samp_n_tot_test,n.samp_n_max))

        storeTreeView.save(taxaTree)


    def selectMicFamilies(self):
        taxaTree = self.taxaTree
        rootNode = taxaTree.getRootNode()
        rootNode.setIndex()
        fgNodes = []
        bgNodes = []
        micNodes = [ taxaTree.getNode(id) for id in micTaxids ]
        for node in taxaTree.iterDepthTop():
            if node.whichSupernode(micNodes) is not None:
                if node.level == self.labelLevel:
                    if node.is_class:
                        if node.samp_n_tot_train >= 80:
                            fgNodes.append(node)
                        else:
                            bgNodes.append(node)
        trainTarg = 30000/len(fgNodes)
        for node in fgNodes:
            node.samp_n_train_targ = min(max(trainTarg,node.sampSel.test.train.min),
                node.samp_n_tot_train)
        bgTrainTarg = trainTarg / 2 / len(bgNodes)
        for node in bgNodes:
            node.samp_n_train_targ = min(bgTrainTarg,node.samp_n_tot_train)
        eukBgNodes = ( taxaTree.getNode(diatomsTaxid), )
        for node in eukBgNodes:
            node.samp_n_train_targ = min(trainTarg/2,node.samp_n_tot_train)
        bgNodes.extend(eukBgNodes)
        return Struct(fgNodes=fgNodes,bgNodes=bgNodes)


    def writeTraining(self,trNodes):
        """Write training data set for each eligible node.
        @pre Nodes have the following attributes: 
        is_class - True for nodes that can serve as classes in SVM
        n_mclass - number of trainable sub-classes
        samp_n_tot_test - number of testing samples in a subtree
        samp_n_tot_train - number of training samples in a subtree (excluding testing)
        is_test - True if a subtree is used for testing
        @post Data sets ready for SVM training are saved. This is done through the
        creation of PredictorDb and descendants of Predictor class.
        @post taxaTree.getRootNode().setIndex() resets internal node index.
        This method must be called after markTraining() sets needed node attributes.
        """
        taxaTree = self.taxaTree
        rootNode = taxaTree.getRootNode()
        rootNode.setIndex()
        self.sampNTrainTot = rootNode.makePropArrayFromAttrib(dtype='i4',name='samp_n_tot_train')
        #for node in taxaTree.iterDepthTop():
        #    assert self.sampNTrainTot[node.lind] == node.samp_n_tot_train
        ## @todo make conditional on not is_test
        def samp_n_choice(node):
            if node.is_test:
                return 0
            else:
                return node.samp_n
        self.sampNTrainSelf = rootNode.makePropArrayFromFunc(dtype='i4',func=samp_n_choice)
        makedir(self.predictorDir)
        kmerReader = KmerReader(hdfSampFile=self.hdfSampFile,sampLen=self.sampLen,kmerLen=self.kmerLen,spacer='N'*self.kmerLen)
        #kmerReader = KmerReaderComb(hdfSampFile=self.hdfSampFile,sampLen=self.sampLen,kmerLenRange=[2,self.kmerLen],spacer='N')
        excludeBranches= [ taxaTree.getNode(id) for id in (35325,35237) ] #dsRNA, dsDNA
        start = time.time()
        if self.kmerRepr == "Sequences":
            svmWriter = SvmStringFeatureWriterTxt(os.path.join(self.predictorDir,self.svmTrainFile))
            svmWriterSplits = (SvmStringFeatureWriterTxt(os.path.join(self.predictorDir,self. svmTrainFile+'.s.1')),SvmStringFeatureWriterTxt(os.path.join(self.predictorDir,self. svmTrainFile+'.s.2')))

        else:
            svmWriter = SvmSparseFeatureWriterTxt(os.path.join(self.predictorDir,self.svmTrainFile))
        #kmerReader.openOutput(predictor.trainFile)
        children = trNodes.fgNodes
        #children = [ n for n in children if sum( (n.isSubnode(e) for e in excludeBranches) ) == 0 ]
        #children = [ n for n in node.iterDepthTop() if n.is_class and \
        #        n.id in ([2,2157]+list(viralTaxidLev2))] #and n.id in (10811,10780,10662,10699)]
        bg_nodes = trNodes.bgNodes
        ids = numpy.zeros(len(children)+2,dtype='i4')
        label = 2
        for child in children:
            child.label = label
            ids[label] = child.id
            label += 1
        for bg_node in bg_nodes:
            bg_node.label = 1
            bg_node.bg_label = 1
        dumpObj(ids,os.path.join(self.predictorDir,"labelToId.pkl"))
        children.extend(bg_nodes)
        for child in children:
            child.splitCounts = numpy.zeros(2,'i4')
        print "Labels-Taxids: ", ids
        # log and pickle the sequences taxids and sample counts for each 'child'
        taxaSampLog = {}
        iChild = 1
        for child in children:
        #for child in node.getChildren(): 
            print "Writing node %s out of %s for label %s" % (iChild,len(children),child.label)
            if child.is_class:
                taxaSampChild = self.pickRandomSampleCounts(child)
                dryRun = False
                if dryRun:
                    nSelected = sum(taxaSampChild.values())
                    if True or nSelected != child.sampSel.test.train.max:
                        print "DEBUG: taxaSampChild = ", \
                                sum(taxaSampChild.values()),  \
                                zip(taxaSampChild.itervalues(), \
                                [ (id,self.sampNTrainTot[node.lind],self.sampNTrainSelf[node.lind]) for (id,node) in 
                                    ( (id,taxaTree.getNode(id)) for id in taxaSampChild.iterkeys() ) ])
                        #pdb.set_trace()
                    continue
                taxaSampLog[child.id] = taxaSampChild
                # we sort taxids by the physical order of sequence db for better performance
                sampNodes = [taxaTree.getNode(taxid) for taxid in taxaSampChild ]
                sampNodes.sort(key=lambda n: n.lnest)
                nWritten = 0
                for taxidSamp in (n.id for n in sampNodes):
                    #nWritten += kmerReader.randomFrequencesWriteSvmSparseTxt(label=child.label,
                    #        taxidSamples=taxidSamp,
                    #        nSamples=taxaSampChild[taxidSamp])
                    #continue
                    kmerMethod = getattr(kmerReader,"random"+self.kmerRepr)
                    for (feature,total) in kmerMethod(taxid=taxidSamp,
                            nSamples=taxaSampChild[taxidSamp]):
                        sampNode = taxaTree.getNode(taxidSamp)
                        assert sampNode is child or sampNode.isSubnode(child)
                        genNode = sampNode.findRankInLineage("genus")
                        if genNode is not None and not genNode.isUnderUnclass:
                            try:
                                splitId = genNode.splitId
                            except AttributeError:
                                scnt = child.splitCounts
                                # among equal minimal counts, add to the randomly picked
                                splitIdMin = scnt.argmin()
                                splitIdsMin = numpy.where(scnt==scnt[splitIdMin])[0]
                                splitId = nrnd.permutation(splitIdsMin)[0]
                                genNode.splitId = splitId
                            child.splitCounts[splitId]+=1
                            svmWriterSplits[splitId].write(child.label,feature)

                        svmWriter.write(child.label,
                                        feature,
                                        taxidSamp)
                        nWritten += 1
                if nWritten < child.samp_n_train_targ:
                    print "Warning: written samples %s - less than requested %s for label %s taxid %s min samples %s" % \
                            (nWritten,child.samp_n_train_targ,child.label,child.id,child.sampSel.test.train.min)
            iChild += 1
        #kmerReader.closeOutput()
        svmWriter.close()
        for svmWriterSplit in svmWriterSplits:
            svmWriterSplit.close()
        dumpObj(taxaSampLog,os.path.join(self.predictorDir,"taxaSamp.train.pkl"))
        print "DEBUG: Done predictor iteration in %s sec" % (time.time() - start)


    def pickRandomSampleCounts(self,node):
        """Randomly pick counts of available samples from taxonomy ids below the node.
        @return dictionary { taxid -> count }.
        Should be called from writeTraining() which sets up sampNTrainTot and sampNTrainSelf
        node property arrays."""
        taxaCounts = {}
        taxaTree = self.taxaTree
        rootNode = taxaTree.getRootNode()
        # Here we rely on 'startIndex' been zero when rootNode.setIndex(startIndex) was called.
        # We could select just a slice self.sampNTrainTot[node.lind:node.rind], but memory wise it 
        # should not matter much because near the tree root the slice would be as big as the
        # entire tree anyway, and indexing is faster with zero based full array.
        # At each call, we create these copies fresh (which is very fast) and keep subtracting
        # from them:
        nTrainTot = self.sampNTrainTot.copy()
        nTrainSelf = self.sampNTrainSelf.copy()
        def actor(n):
            choices = [ child for child in n.children if nTrainTot[child.lind] > 0 ]
            choicesLongSeq = [ child for child in choices if child.has_long_seq ]
            if len(choicesLongSeq) > 0:
                choices = choicesLongSeq
            nChoices = len(choices)
            lind = n.lind
            if nTrainSelf[lind] > 0 and (n.has_long_seq_self or len(choicesLongSeq) == 0):
                # random_integers returns closed range [low,high]
                iChoice = nrnd.random_integers(0,nChoices) - 1 # -1 <= iChoice < len(choices)
                if iChoice < 0:
                    nTrainSelf[lind] -= 1
                    nTrainTot[lind] -= 1
                    try:
                        taxaCounts[n.id] += 1
                    except KeyError:
                        taxaCounts[n.id] = 1
                    return
                # If we got here, iChoice points to one of the subnodes
            elif nChoices > 0:
                iChoice = nrnd.randint(0,nChoices)
            else:
                raise ValueError("Tree error")
            actor(choices[iChoice])
            nTrainTot[lind] -= 1

        nTrainTarget = node.samp_n_train_targ
        nTrainDone = 0
        while nTrainDone < nTrainTarget:
            assert nTrainTot[node.lind] > 0
            actor(node)
            nTrainDone += 1
        return taxaCounts


    def writeTesting(self):
        """Write testing data set.
        Currently, we write one dataset and use it on all classifiers.
        The reason is because we use external executables for prediction,
        and starting the executable for different models along the tree
        path for each sample would be too expensive. Instead, we subsample
        all testing taxa and create a single file. Along the way, we
        also create an HDF dataset with a lineage of predictors for 
        each unique taxid among selected samples. We also save the taxid
        of each sample in an HDF array.
        @pre Nodes have the following attributes: 
        is_class - True for nodes that can serve as classes in SVM
        n_mclass - number of trainable sub-classes
        samp_n_tot_test - number of testing samples in a subtree
        is_test - True if a subtree is used for testing
        This method must be called after markTraining() sets needed node attributes.
        """
        taxaTree = self.taxaTree
        rootNode = taxaTree.getRootNode()
        rootNode.setIndex()
        class SampTest(pt.IsDescription):
            lind = pt.Int32Col(pos=1)
        hdfFile = pt.openFile(os.path.join(self.predictorDir,self.hdfTestFile), mode='w')
        hdfSamples = hdfFile.createTable(hdfFile.root, 'sampTest', SampTest, "Testing Samples")
        sampTestDtype = hdfDtype(hdfSamples) 
        kmerReader = KmerReader(hdfSampFile=self.hdfSampFile,
                sampLen=self.sampLen,
                kmerLen=self.kmerLen,
                spacer='N')
        #kmerReader = KmerReaderComb(hdfSampFile=self.hdfSampFile,
        #        sampLen=self.sampLen,
        #        kmerLenRange=[2,self.kmerLen],
        #        spacer='N')
        kmerReader.setSubSamplerUniRandomEnd(minLen=self.minTestSampLen,maxLen=self.maxTestSampLen)
        #kmerReader.openOutput(self.getFileName(self.svmTestFile))
        if self.kmerRepr == "Sequences":
            svmWriter = SvmStringFeatureWriterTxt(os.path.join(self.predictorDir,self.svmTestFile))
        else:
            svmWriter = SvmSparseFeatureWriterTxt(os.path.join(self.predictorDir,self.svmTestFile))
        taxaSampLog = {}
        children = [ n for n in self.taxaTree.iterDepthTop() if hasattr(n,"label") ]
        for child in children:
            taxaSampLog[child.id] = {}
            for node in child.iterDepthTop():
                if node.is_test and node.samp_n > 0:
                    nSampTestTarg = min(node.samp_n,self.maxTestSampPerTaxa)
                    #DEBUG TMP:
                    #lin = node.lineage()
                    #print [ n.name for n in reversed(lin) ]
                    #nSampTestOut = kmerReader.randomFrequencesWriteSvmSparseTxt(label=child.label,
                    #        taxidSamples=node.id,
                    #        nSamples=nSampTestTarg)
                    nSampTestOut = 0
                    kmerMethod = getattr(kmerReader,"random"+self.kmerRepr)
                    for (feature,total) in kmerMethod(taxid=node.id,
                                nSamples=nSampTestTarg):
                        assert node == child or node.isSubnode(child)
                        svmWriter.write(child.label,
                                        feature,
                                        node.id)
                        nSampTestOut += 1
                    samples = numpy.zeros(nSampTestOut,dtype=sampTestDtype)
                    samples['lind'] = node.lind
                    hdfSamples.append(samples)
                    taxaSampLog[child.id][node.id] = nSampTestOut
        #kmerReader.closeOutput()
        kmerReader.clearSubSampler()
        svmWriter.close()
        dumpObj(taxaSampLog,os.path.join(self.predictorDir,"taxaSamp.test.pkl"))
        hdfSamples.flush()

    def tmp_hdfSamplesConv(self):

        taxaTree = self.taxaTree
        rootNode = taxaTree.getRootNode()
        rootNode.setIndex()
        hdfFile = pt.openFile(os.path.join(self.predictorDir,self.hdfTestFile), mode='r')
        hdfSamples = hdfFile.getNode(hdfFile.root, 'sampTest')
        sampTestDtype = hdfDtype(hdfSamples)
        sampLind = hdfSamples.read(field='lind')
        lindToId = rootNode.makePropArrayFromAttrib(dtype='i4',name='id')
        sampId = lindToId[sampLind]
        dumpObj(sampId,os.path.join(self.predictorDir,"sampTaxid.test.pkl"))
        # make genus index and mark 2,3,4 splits as columns in Nx3 array where N is sample count
        cnt = numpy.bincount(sampId)
        idsGen = numpy.zeros_like(cnt)
        for id in numpy.where(cnt>0)[0]:
            nodeGen = taxaTree.getNode(id).findRankInLineage(rank="genus")
            if nodeGen is not None:
                idsGen[id] = nodeGen.id
        sampIdGen = idsGen[sampId]
        gen = numpy.unique(sampIdGen)
        gen = gen.compress(gen>0)
        nSplitSteps = 3 # will be splitted into (2,3,4) parts
        sampSplit = numpy.zeros((len(sampIdGen),nSplitSteps),dtype='i2')
        for iSplitStep in xrange(nSplitSteps):
            nSplit = 2 + iSplitStep
            genSplit = nrnd.random_integers(1,nSplit,len(gen))
            sampSplit[:,iSplitStep] = genSplit[gen.searchsorted(sampIdGen)]
            sampSplit[sampIdGen == 0,iSplitStep] = 0
        dumpObj(sampSplit,os.path.join(self.predictorDir,"sampSplit.test.pkl"))
        #pdb.set_trace()

    def writeTestingTree(self):
        """TMP: Write HDF file with tree structure data needed for loading testing results."""
        taxaTree = self.taxaTree
        # (id,rind)[lind]

    def submitTraining(self):
        "TMP: Submit training batch jobs"

        btOpt = self.batchRun

        bt = BatchSubmitter(PROJECT_CODE=btOpt.PROJECT_CODE,MEM=btOpt.MEM,ARCH=btOpt.ARCH)

        iPredictor = 0
        for predDir in [ os.path.join(self.predictorDir,d) for d in sorted(os.listdir(self.predictorDir)) ]:
            svm = SVMLibLinear(workDir=os.path.join(predDir,"ms"))
            bt.submitIf(maxQueued=btOpt.maxQueued,cmd=' '.join(svm.trainCmd(trainFile=os.path.join(predDir,"train.svm"))),
                    scriptName="train",runDir=svm.workDir,sleepTime=1,dryRun=False)
            print "DEBUG: iPredictor = ", iPredictor
            iPredictor += 1
            if iPredictor > 0:
                break

    
    def mkDbTestingTaxaDb(self,rank="genus",superRank="family"):
        """Select complete 'rank' nodes uniformly distributed across 'superRank' super-nodes."""
        rankId = self.levels.getLevelId(rank)
        superRankId = self.levels.getLevelId(superRank)
        tblTxSamp = self.getTableName("tx_samp")
        """
        select id 
        from %(tblTree)s tr_leaf, %(tblTxSamp) tx_samp, %(tblTree)s tr_sup
        where 
        tr_leaf.id = tx_samp.id and
        %(subtreeLeafSup)s and
        tr_sup.idlevel >= %(superRankId)i and min_depth
        """ % { 'subtreeLeafSup' : sqlIsSubTree('tr_leaf', 'tr_sup'),
                'superRankId'    : superRankId }
            

class RandomSampPicker:
    pass
