### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Applicaion for  building and using IMM-based taxonomic classifier in the spirit of Phymm"""

from MGT.ImmApp import *
from MGT.Taxa import *
from MGT.App import *
from MGT.DirStore import *
from MGT.SeqDbFasta import *
from MGT.ArchiveApp import *

class TaxaPred(Struct):
    """Class to represent taxonomic predictions - one taxid per sample"""
    pass

class ImmClassifierApp(App):
    """App-derived class for building and using IMM-based taxonomic classifier in the spirit of Phymm.
    The main difference with Phymm is that we build IMMs for higher-level clades by pulling sequence
    data for subnodes.
    This class can be mostly viewed as imposing a TaxaTree structure onto ImmApp."""

    batchDepModes = ("predict","score","train")

    ## Special taxonomy ID value to mark rejected samples 
    rejTaxid = 0

    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        
        optChoicesMode = ("score","predict","train","setup-train","proc-scores")
        
        def optParseCallback_StoreAbsPath(option, opt_str, value, parser):
            setattr(parser.values, option.dest, os.path.abspath(value))

        def optParseMakeOption_Path(shortName,longName,dest,help=None,default=None):
            return make_option(shortName,longName,
            action="callback", 
            callback=optParseCallback_StoreAbsPath,
            type="string",
            dest=dest,
            default=default,
            help=help)
        
        option_list = [
            
            make_option("-m", "--mode",
            action="store", 
            type="choice",
            choices=optChoicesMode,
            dest="mode",
            default="predict",
            help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            
            optParseMakeOption_Path(None, "--db-seq",
            dest="seqDb",
            help="SeqDbFasta path, used either to create custom SeqDb from scaffolds "+\
                "or to use an existing SeqDb for training the IMMs"),
            
            make_option(None, "--db-imm",
            action="append", 
            type="string",
            dest="immDb",
            help="Path to a collection of IMMs stored as a directory. "+\
                    "Multiple entries are allowed in prediction mode, "+\
                    "but only one entry - in training mode."),
            
            make_option(None, "--db-imm-archive",
            action="append", 
            type="string",
            dest="immDbArchive",
            help="Similar to --db-imm, but each IMM collection is represented by a single"+\
                    "archive file. This can be used both in training and prediction phases. "+\
                    "Multiple entries are allowed in prediction mode. "+\
                    "The list of collections defined with this option will be concatenated "+\
                    "with the list defined with --db-imm option."),
            
            optParseMakeOption_Path(None, "--imm-seq-ids",
            default="imm-seq-ids",
            dest="immIdToSeqIds",
            help="File that maps IMM IDs to lists of seq IDs during IMM training"),
            
            optParseMakeOption_Path(None, "--imm-ids",
            dest="immIds",
            help="File with list of IMM IDs to use in scoring and prediction. Default is all"+\
                    " IMMs from --imm-seq-ids"),
            
            optParseMakeOption_Path(None, "--inp-seq",
            dest="inpSeq",
            help="File with input FASTA sequence for prediction"),

            optParseMakeOption_Path(None, "--inp-train-seq",
            dest="inpTrainSeq",
            help="File with input FASTA sequences for training the models"),
            
            make_option(None, "--max-seq-id-cnt",
            action="store", 
            type="int",
            default=100,
            dest="maxSeqIdCnt",
            help="Maximum number of training SeqDB IDs to propagate up from "+\
                    "every child of a given node"),
            
            make_option(None, "--n-imm-batches",
            action="store", 
            type="int",
            default=200,
            dest="nImmBatches",
            help="Try to split processing into that many batches for each ICM set "+\
                    "(leading to separate jobs in batch run-mode)"),
            
            optParseMakeOption_Path(None, "--out-dir",
            default="results",
            dest="outDir",
            help="Directory name for output score files"),
            
            optParseMakeOption_Path(None, "--out-score-comb",
            dest="outScoreComb",
            help="Output file for combined raw scores. Default "+\
                    "is 'combined'+ImmApp.scoreSfx inside --out-dir"),
        
            optParseMakeOption_Path(None, "--taxa-tree-pkl",
            dest="taxaTreePkl",
            help="Custom taxonomy tree saved in pickle format. If not set, standard NCBI tree is used, "+\
                    "except when training custom IMMs when this has to be always defined."),
            
            make_option(None, "--pred-min-len-samp",
            action="store", 
            type="int",
            dest="predMinLenSamp",
            help="Min length of samples to consider for prediction. 300 is default "+\
                    "for bacterial classification; 5000 is default for predicting "+\
                    "hosts for viral contigs."),
            
            make_option(None, "--train-min-len-samp",
            action="store", 
            type="int",
            default=100000,
            dest="trainMinLenSamp",
            help="Min length of scaffolds to consider for training custom IMMs"),
            
            make_option(None, "--new-tax-name-top",
            action="store", 
            type="string",
            default="mgt_reference",
            dest="newTaxNameTop",
            help="Root name for custom reference sequences"),
            
            optParseMakeOption_Path(None, "--pred-out-dir",
            default="results",
            dest="predOutDir",
            help="Output directory for classification results"),
            
            optParseMakeOption_Path(None, "--pred-out-taxa",
            dest="predOutTaxa",
            help="Output file with predicted taxa; default is pred-taxa inside "+\
                    "--pred-out-dir"),
            
            optParseMakeOption_Path(None, "--pred-out-taxa-csv",
            dest="predOutTaxaCsv",
            help="Output CSV file with predicted taxa; default is --pred-out-taxa "+\
                    "+ '.csv'"),
            
            make_option(None, "--rej-ranks-higher",
            action="store", 
            type="string",
            dest="rejectRanksHigher",
            help="If a sample was assigned a clade above this rank or below , it will be marked as "+\
                    "'rejected' instead. The default value of 'superkingdom' effectively disables "+\
                    "this filer. Set to 'order' if performing phage-host assignment."),
            
            make_option(None, "--pred-mode",
            action="store",
            type="choice",
            choices=("host","taxa"),
            default="taxa",
            dest="predMode",
            help="Set the prediction mode: 'host' will work in a mode that assigns "+\
                    "host taxonomy to (presumed) bacteriophage "+\
                    "sequence. Setting this will overwrite the value of some other "+\
                    "options (currently it will set --rej-ranks-higher=order and "+\
                    "--pred-min-len-samp=5000 if they were not defined. "+\
                    "'bac' [default] will try to assign bacterial taxonomy to the "+\
                    "input sequences."),
        ]
        return Struct(usage = klass.__doc__+"\n"+\
                "%prog [options]",option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.setIfUndef("cwd",os.getcwd())
        if ( not opt.immDbArchive and not opt.immDb ):
            opt.immDb = [ pjoin(opt.cwd,"imm") ]
        opt.setIfUndef("seqDb",pjoin(opt.cwd,"seqDb"))
        opt.setIfUndef("immIds",opt.immIdToSeqIds)
        opt.setIfUndef("outScoreComb",pjoin(opt.outDir,"combined"+ImmApp.scoreSfx))
        opt.setIfUndef("predOutTaxa",pjoin(opt.predOutDir,"pred-taxa"))
        opt.setIfUndef("predOutTaxaCsv",opt.predOutTaxa+".csv")
        opt.setIfUndef("predOutDbSqlite",opt.predOutTaxa+".sqlite")
        opt.setIfUndef("predOutStatsDir", pjoin(opt.predOutDir,"stats"))
        opt.setIfUndef("predOutStatsCsv", pjoin(opt.predOutStatsDir,"stats.csv"))
        opt.setIfUndef("newTaxidTop",mgtTaxidFirst)
        opt.setIfUndef("immDbWorkDir",pjoin(opt.cwd,"immDbWorkDir"))
        opt.setIfUndef("scoreWorkDir",pjoin(opt.cwd,"scoreWorkDir"))
        if opt.predMode == "host": 
            opt.setIfUndef("rejectRanksHigher","order")
            opt.setIfUndef("predMinLenSamp",5000)
        elif opt.predMode == "taxa":
            opt.setIfUndef("rejectRanksHigher","superkingdom")
            opt.setIfUndef("predMinLenSamp",300)
        else:
            raise ValueError("Unknown --pred-mode value: %s" % (opt.predMode,))

    
    def instanceOptionsPost(self,opt):
        """Set (in place) instance-specific options.
        This is called from __init__() and has access to the execution context (such as current dir)."""
        ## parseCmdLinePost will not modify options that are already defined, so we need to do it here
        if isinstance(opt.immDb,str):
            opt.immDb = [ opt.immDb ]
        elif opt.immDb is None:
            opt.immDb = list()
        if isinstance(opt.immDbArchive,str):
            opt.immDbArchive = [ opt.immDbArchive ]
        elif opt.immDbArchive is None:
            opt.immDbArchive = list()
    
    def initWork(self,**kw):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.taxaLevels = None #will be lazy-loaded
        self.seqDb = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        makedir(opt.immDbWorkDir)
        makedir(opt.scoreWorkDir)
   

    def doWork(self,**kw):
        """@todo We currently cannot run train and score in one batch
        submission because score stage needs train results to figure
        out the IMM list. It will be easy to fix by saving list of
        IMM IDs during train submission and reusing it during score
        submission. A more economical solution would be the CA way -
        submitting a dumb terminator job, submitting the correct 
        scoring jobs at the end of training job, and using qmod to
        make terminator job dependent on those new scoring jobs.
        This approach would allow making initial user submission 
        process very light and fast, while currently a fairly heavy
        expansion and reading of training FASTA file is required
        during submission for training custom IMMs. The downside
        is an extra requirement on the LRM config that it has to
        allow job submission from compute nodes."""
        opt = self.opt
        ret = None
        if opt.mode == "train":
            ret = self.train(**kw)
        elif opt.mode == "score":
            ret = self.score(**kw)
        elif opt.mode == "predict":
            ret = self.predict(**kw)
        elif opt.mode == "proc-scores":
            ret = self.processImmScores(**kw)
        elif opt.mode == "setup-train":
            return self.setupTraining(**kw)
        elif opt.mode == "combine-scores":
            ret = self.combineScores(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))
        return ret

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree(pklFile=self.opt.taxaTreePkl)
        return self.taxaTree

    def getTaxaLevels(self):
        if self.taxaLevels is None:
            #that assigns "level" and "idlevel" attributes to TaxaTree nodes
            self.taxaLevels = TaxaLevels(self.getTaxaTree())
        return self.taxaLevels
    
    def getSeqDb(self):
        opt = self.opt
        if self.seqDb is None:
            self.seqDb = SeqDbFasta.open(opt.seqDb,mode="r") #"r"
            return self.seqDb


    ## Methods that assign training sequences to higher-level nodes

    def mapSeqToTree(self):
        """Assign list of SeqDB IDs to corresponding taxonomy tree nodes.
        In the current SeqDB version, each ID is a unique taxonomy id, so
        a one-element list will be assigned.
        @post Attribute 'leafSeqDbIds' is assigned to EVERY node and contains a list of IDs, possibly empty.
        The empty list will be a reference to a single shared object, so it should be treated as immutable"""
        taxaTree = self.getTaxaTree()
        seqDb = self.getSeqDb()
        taxaList = seqDb.getTaxaList()
        emptyList = []
        taxaTree.setAttribute("leafSeqDbIds",emptyList,doCopy=False)
        for taxid in taxaList:
            taxaTree.getNode(taxid).leafSeqDbIds = [ taxid ]


    def pickSeqOnTree(self,maxSeqIdCnt):
        """Assign to each node of the taxonomy tree the list of SeqDB IDs in the subtree under that node.
        It picks a medium count of IDs from each child node and itself, clipped by maxSeqIdCnt.
        """
        taxaTree = self.getTaxaTree()

        def actor(node):
            if node.isLeaf() and len(node.leafSeqDbIds)>0:
                node.pickedSeqDbIds = node.leafSeqDbIds
            else:
                chSeqIds = [c.pickedSeqDbIds for c in node.getChildren() if hasattr(c,"pickedSeqDbIds")]
                chSeqIds.append(node.leafSeqDbIds)
                chSeqIds = [ x for x in chSeqIds if len(x) > 0 ]
                if len(chSeqIds) > 0:
                    targLen = int(n.median([ len(x) for x in chSeqIds ]))
                    targLen = min(maxSeqIdCnt,targLen)
                    chSeqIds = [ sampleBoundWOR(x,targLen) for x in chSeqIds ]
                    node.pickedSeqDbIds = sum(chSeqIds,[])
                    assert len(node.pickedSeqDbIds) > 0

        taxaTree.visitDepthBottom(actor)

    def defineImms(self):
        taxaTree = self.getTaxaTree()
        immIdToSeqIds = {}
        for node in taxaTree.iterDepthTop():
            if hasattr(node,"pickedSeqDbIds"):
                immIdToSeqIds[node.id] = node.pickedSeqDbIds
        dumpObj(immIdToSeqIds,self.opt.immIdToSeqIds)

    def makeCustomTaxaTreeAndSeqDb(self,**kw):
        """Add provided scaffolds as custom nodes to the taxonomy tree and save the tree, and also create SeqDbFasta for them.
        The resulting SeqDbFasta and TaxaTree can be used to train the IMM models.
        Each scaffold is treated as a separate taxonomic unit under the super-node of environmental bacterial sequences."""
        opt = self.opt
        assert opt.taxaTreePkl, "File name for a custom taxonomic tree must be provided"
        self.seqDb = None
        seqDb = SeqDbFasta.open(path=opt.seqDb,mode="c")
        compr = SymbolRunsCompressor(sym="N",minLen=1)
        nonDegenSymb = "ATCGatcg"
        newTaxaTop = TaxaNode(id=opt.newTaxidTop,name=opt.newTaxNameTop,
                rank=unclassRank,divid=dividEnv,names=list())
        nextNewTaxid = newTaxaTop.id + 1
        fastaReader = FastaReader(opt.inpTrainSeq)
        nNodesOut = 0
        for rec in fastaReader.records():
            hdr = rec.header()
            seqid = rec.getId() # scf768709870
            seq = compr(rec.sequence())
            if not checkSaneAlphaHist(seq,nonDegenSymb,minNonDegenRatio=0.99):
                print "WARNING: ratio of degenerate symbols is too high, "+\
                        "skipping the reference scaffold id %s" % (seqid,)
            if len(seq) >= opt.trainMinLenSamp:
                taxaNode = TaxaNode(id=nextNewTaxid,name=opt.newTaxNameTop+"_"+seqid,rank=unclassRank,divid=dividEnv,names=list())
                taxaNode.setParent(newTaxaTop)
                # that will be used by the following call to defineImms()
                taxaNode.pickedSeqDbIds = [ taxaNode.id ]
                fastaWriter = seqDb.fastaWriter(id=nextNewTaxid,lineLen=80)
                fastaWriter.record(header=hdr,sequence=seq)
                fastaWriter.close()
                nextNewTaxid += 1
                nNodesOut += 1
        print "DEBUG: written %s sequence files in SeqDbFasta" % (nNodesOut,)
        if nNodesOut <= 0:
            raise ValueError(("No training nodes conform to the minimum sequence length "+\
                    "requirement of %s (after compressing degenerate runs)") % (opt.trainMinLenSamp,))
        fastaReader.close()
        self.taxaTree = None
        taxaTree = loadTaxaTree() # pristine NCBI tree
        envBacTop = taxaTree.getNode(envBacTaxid)
        newTaxaTop.setParent(envBacTop)
        taxaTree.rebuild()
        taxaTreeStore = NodeStoragePickle(opt.taxaTreePkl)
        taxaTreeStore.save(taxaTree)
        self.taxaTree = taxaTree
    
    def setupTraining(self,**kw):
        opt = self.opt
        if opt.inpTrainSeq:
            #this also sets node.pickedSeqDbIds,
            #and only those nodes will be used to train
            #models, w/o propagation up the tree, because
            #otherwise we would get intersecting set of models
            #with the standard reference set
            self.makeCustomTaxaTreeAndSeqDb()
        else:
            self.mapSeqToTree()
            self.pickSeqOnTree(opt.maxSeqIdCnt)
        self.defineImms()

    def _archiveNameToDirName(self,archiveName,topDir,subDir=None):
        import tempfile
        d = tempfile.mkdtemp(suffix=".tmp",
                prefix=os.path.basename(archiveName),
                dir=topDir)
        if subDir is not None:
            d = pjoin(d,subDir)
            makedir(d)
        return d

    def _archiveNamesToDirNames(self,archiveNames,topDir,subDir=None):
        return [ (self._archiveNameToDirName(arch,topDir,subDir),arch) for \
                arch in archiveNames ]
    
    def _immDbNameToScoreDirName(self,immDbName,topDir):
        import tempfile
        return tempfile.mkdtemp(suffix=".tmp",
                prefix=os.path.basename(immDbName),
                dir=topDir)
    
    def train(self,**kw):
        """Train all IMMs.
        Parameters are taken from self.opt
        """
        opt = self.opt

        immDb = [ (d,None) for d in opt.immDb ]
        immDbArch = self._archiveNamesToDirNames(opt.immDbArchive,opt.immDbWorkDir)
        immDb += immDbArch
        #changing this will also require separate setupTraining() for each SeqDb etc
        assert len(immDb) == 1,"Only one IMM DB is allowed during training: %s" % (immDb,)

        optI = copy(opt)
        optI.mode = "setup-train"
        #The SeqDb must be available before training is submitted because
        #the IDs of training models are decided during submission.
        #TODO: figure out some way to change the situation described above,
        #although it might be difficult. As of now, running setup-train 
        #inproc can be lengthy.
        optI.runMode = "inproc"
        app = self.factory(opt=optI)
        jobs = app.run(**kw)
        
        optI = copy(opt)
        optI.mode = "train"

        optI.immDb = immDb[0][0]

        app = ImmApp(opt=optI)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = app.run(**kw)

        if immDbArch:
            dbArch = immDb[0]
            optI = copy(opt)
            optI.mode = "archive"
            optI.path = pjoin(dbArch[0],"") #append '/' to make tar-bomb
            optI.archive = dbArch[1]
            optI.safe = True # currently noop in mode="archive"
            app = ArchiveApp(opt=optI)
            kw = kw.copy()
            kw["depend"] = jobs
            jobs = app.run(**kw)
        
        return jobs

    def score(self,**kw):
        """Score with all IMMs.
        Parameters are taken from self.opt
        @param inpSeq Name of the input multi-FASTA file to score
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        immDb = [ (d,None) for d in opt.immDb ]
        immDbArch = self._archiveNamesToDirNames(opt.immDbArchive,opt.immDbWorkDir,"imm")
        immDb += immDbArch
        jobsD = kw.get("depend",list())
        jobs = []
        outSubScores = []
        for immD in immDb:
            jobsI = copy(jobsD)
            if immD[1] is not None:
                optI = copy(opt)
                optI.mode = "extract"
                optI.path = immD[0]
                optI.archive = immD[1]
                optI.safe = True
                app = ArchiveApp(opt=optI)
                #optI.path is not actually used below
                immIds = ImmStore(immD[0]).listImmIds(iterPaths=(item.name for item in app.iterMembers()))
                kwI = kw.copy()
                jobsI = app.run(**kwI)
            else:
                immIds = ImmStore(immD[0]).listImmIds()
            assert len(immIds) > 0,"No IMMs found in IMM DB - probably training did not run yet"
            
            kwI = kw.copy()
            kwI["depend"] = jobsI
        
            optI = copy(opt)
            optI.mode = "score"
            optI.immDb = immD[0]
            #Until we change ImmApp, we have to run each in a separate dir
            #because it uses a fixed file name for immIds file that it generates.
            optI.cwd = self._immDbNameToScoreDirName(optI.immDb,opt.scoreWorkDir)
            optI.outDir = optI.cwd
            #TODO: have a separate set of the options below for each immDb,
            #until then, we have to zap them and cause ImmApp to use all
            #available imms in each collection.
            optI.immIdToSeqIds = None
            optI.immIds = pjoin(optI.outDir,"imm-ids.pkl")
            dumpObj(immIds,optI.immIds)
            optI.outScoreComb = pjoin(optI.outDir,"combined"+ImmApp.scoreSfx)
            outSubScores.append(optI.outScoreComb)
            app = ImmApp(opt=optI)
            jobsI = app.run(**kwI)
            jobs += jobsI
        
        optI = copy(opt)
        optI.mode = "combine-scores"
        optI.outSubScores = outSubScores
        
        app = self.factory(opt=optI)
        kwI = kw.copy()
        kwI["depend"] = jobs
        jobs = app.run(**kwI)
        
        return jobs

    def combineScores(self,**kw):
        """Combine scores from different immDbs.
        Parameters are taken from self.opt
        @param outSubScores List of files with input score matrices
        @param outScoreComb Name for output file with combined scores
        """
        opt = self.opt
        subScores = [ loadObj(subScoreFile) for subScoreFile in opt.outSubScores ]
        # calling class method:
        scoreComb = subScores[0].catImms(subScores)
        makeFilePath(opt.outScoreComb)
        dumpObj(scoreComb,opt.outScoreComb)

    def _maskScoresNonSubtrees(self,taxaTree,immScores,posRoots):
        """Set to a negative infinity (numpy.NINF) all columns in score matrix that point to NOT subtrees of posRoots nodes.
        This is used to mask all scores pointing to other than bacteria or archaea."""
        scores = immScores.scores
        idImms = immScores.idImm
        for (iCol,idImm) in enumerate(idImms):
            #assume here that idImm is a taxid
            node = taxaTree.getNode(idImm)
            if not (node.isSubnodeAny(posRoots) or node in posRoots):
                scores[:,iCol] = n.NINF

    def _maskScoresByRanks(self,taxaTree,immScores):
        """Set to a negative infinity (numpy.NINF) all columns in score matrix that point to nodes of ranks not in desired range.
        """
        opt = self.opt
        scores = immScores.scores
        idImms = immScores.idImm
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        max_linn_levid = taxaLevels.getLevelId(opt.rejectRanksHigher)
        #min_linn_levid = max_linn_levid
        # Safe to set it to 0, because there is also at least superkingdom above 
        # min_linn_levid = taxaLevels.getLinnLevelIdRange()[0]
        min_linn_levid = 0
        # this still can exclude nodes that are above subtree that has no IMMs but 
        # below lowest max_lonn_levid (e.g. ref sequence attached directly to sub-species
        # that has strain nodes w/o sequence, and family node immediately above.
        # We need to assign is_leaf_imm attribute to taxaTree nodes to do it right.
        for (iCol,idImm) in enumerate(idImms):
            #assume here that idImm is a taxid
            node = taxaTree.getNode(idImm)
            # node.isLeaf() takes care of env sequences and ref strains under family
            if  (not node.isLeaf()) and \
                (not taxaLevels.isNodeInLinnLevelRange(node,min_linn_levid,max_linn_levid)):
                scores[:,iCol] = n.NINF
   
    def _taxaTopScores(self,taxaTree,immScores,topScoreN):
        """Get TaxaTree nodes of topScoreN scores for each sample.
        """
        scores = immScores.scores
        idImms = immScores.idImm
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        #get indices of topScoreN largest scores in each row
        indScores = scores.argsort(axis=1)[:,-1:-topScoreN-1:-1]
        return idImms[indScores]

    def _getImmScores(self,reload=False):
        if not hasattr(self,"immScores") or reload:
            self.immScores = loadObj(self.opt.outScoreComb)
            self.immScores.idImm = n.asarray(self.immScores.idImm,dtype=int)
        return self.immScores


    def _normalizeScores(self,immScores,immScoresRnd):
        sc = immScores
        scRand = immScoresRnd
        scoreRand = scRand.scores
        baseScoreRand = (scoreRand.T/scRand.lenSamp).T
        meanRand = baseScoreRand.mean(0)
        stdRand = baseScoreRand.std(0)
        normScoreRand = (baseScoreRand - meanRand)/stdRand
        baseScores = (sc.scores.T/sc.lenSamp).T
        assert n.all(scRand.lenSamp==scRand.lenSamp[0])
        ratLen = sc.lenSamp.astype(float)/scRand.lenSamp[0]
        # base score is an average of N random variables,
        # and therefore has a variance var(1-base)/N
        normScore = (baseScores - meanRand)/stdRand
        normScore = (normScore.T * ratLen**0.5).T
        #normScore = baseScores/meanRand
        #normScore = baseScores
        #pdb.set_trace()
        # we could also return dot(stdRand,ratLen**0.5) to give
        # avg base score std per given sample length by a given IMM on random
        # sequence, but the distribution on our random set is virtually normal,
        # and it is much to the left from the distribution of actual samples,
        # so there is no sense in using percentile cutoffs like 95%, which is
        # at about 1.645 z-score for standard normal distribution
        return normScore


    def processImmScores(self,**kw):
        """Process raw IMM scores to predict taxonomy.
        This handles both classification of bacterial sequences and host assignment 
        for viral sequences.
        Parameters are taken from self.opt.
        @param outScoreComb File with ImmScores object
        @param rndScoreComb File with ImmScores object for random query sequences
        @param predOutTaxa Output file with predicted taxa
        """
        opt = self.opt
        sc = loadObj(opt.outScoreComb)
        #assume idImm are str(taxids):
        sc.idImm = n.asarray(sc.idImm,dtype=int)
        #scRnd = loadObj(opt.rndScoreComb)
        #sc.scores = self._normalizeScores(sc,scRnd)
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        scores = sc.scores
        idImms = sc.idImm
        micRoots = tuple([ taxaTree.getNode(taxid) for taxid in micTaxids ])
        virRoot = taxaTree.getNode(virTaxid)
        #scVirRoot = scores[:,idImms == virTaxid][:,0]
        #scCellRoot = scores[:,idImms == cellTaxid][:,0]
        self._maskScoresNonSubtrees(taxaTree,immScores=sc,posRoots=micRoots)
        if opt.rejectRanksHigher is not "superkingdom":
            self._maskScoresByRanks(taxaTree,immScores=sc)
        topTaxids = self._taxaTopScores(taxaTree,immScores=sc,topScoreN=10)
        argmaxSc = scores.argmax(1)
        maxSc = scores.max(1)
        ## @todo in case of score ties, pick the lowest rank if it is a single entry
        predTaxids = idImms[argmaxSc]
        #predTaxids[maxSc<10] = self.rejTaxid
        # this will reject any sample that has top level viral score more
        # that top level cellular org score, on the assumption that easily
        # assignbale viruses should be closer to cellular orgs than to other
        # viruses.
        # Result: removed 90 out of 250 samples, with no change in specificity.
        #predTaxids[scVirRoot>scCellRoot] = self.rejTaxid
        # this excluded 10 out of 430, no change in specificity
        #predTaxids[scVirRoot>=maxSc] = self.rejTaxid
        
        # This not the same as _maskScoresByRanks() above (because this
        # will actually assign a reject label.
        if opt.rejectRanksHigher is not "superkingdom":
            max_linn_levid = taxaLevels.getLevelId(opt.rejectRanksHigher)
            min_linn_levid = taxaLevels.getLinnLevelIdRange()[0]
            # Reject predictions to clades outside of certain clade level range,
            # as well as to any viral node
            # This rejected 36 out of 450 and resulted in 2% improvement in specificity
            for i in xrange(len(predTaxids)):
                if not predTaxids[i] == self.rejTaxid:
                    predNode = taxaTree.getNode(predTaxids[i])
                    if predNode.isUnder(virRoot):
                        predTaxids[i] = self.rejTaxid
                    # we need to protect leaf nodes because we place environmental scaffolds
                    # as no_rank under bacteria->environmental
                    elif not (predNode.isLeaf() or taxaLevels.isNodeInLinnLevelRange(predNode,
                            min_linn_levid,max_linn_levid)):
                        predTaxids[i] = self.rejTaxid

        pred = TaxaPred(idSamp=sc.idSamp,predTaxid=predTaxids,topTaxid=topTaxids,lenSamp=sc.lenSamp)
        makeFilePath(opt.predOutTaxa)
        dumpObj(pred,opt.predOutTaxa)
        self._exportPredictions(taxaPred=pred)
        self._statsPred()
        #todo: "root" tax level; matrix of predhost idlevel's; batch imm.score() jobs; score hist;

    def _exportPredictions(self,taxaPred):
        opt = self.opt    
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        levNames = taxaLevels.getLevelNames("ascend")
        makeFilePath(opt.predOutTaxaCsv)
        out = openCompressed(opt.predOutTaxaCsv,"w")
        flds = ["id","len","taxid","name","rank"]
        for lev in levNames:
            flds += ["taxid_"+lev,"name_"+lev]
        w = csv.DictWriter(out, fieldnames=flds, restval='null',dialect='excel-tab')
        w.writerow(dict([(fld,fld) for fld in flds]))
        for idS,taxid,lenS in it.izip(taxaPred.idSamp,taxaPred.predTaxid,taxaPred.lenSamp):
            if lenS < opt.predMinLenSamp:
                continue
            row = dict(id=idS,len=lenS,taxid=taxid)
            if taxid != self.rejTaxid:
                node = taxaTree.getNode(taxid)
                row["name"] = node.name
                row["rank"] = node.linn_level
                lin = taxaLevels.lineage(node,withUnclass=False)
                for ln in lin:
                    row["name_"+ln.level] = ln.name
                    row["taxid_"+ln.level] = ln.id
            w.writerow(row)
        out.close()
        
        makeFilePath(opt.predOutDbSqlite)
        db = DbSqlLite(dbpath=opt.predOutDbSqlite)
        db.createTableFromCsv(name="scaff_pred",
                csvFile=opt.predOutTaxaCsv,
                hasHeader=True,
                indices={"names":("id","taxid","name",)})
        db.close()

    def _statsPred(self,**kw):
        """Create aggregate tables and csv files to show various relationships between predictions"""
        opt = self.opt
        db = DbSqlLite(dbpath=opt.predOutDbSqlite)
        rmrf(opt.predOutStatsDir)
        makeFilePath(opt.predOutStatsCsv)
        outCsv = openCompressed(opt.predOutStatsCsv,'w')
        sqlAsComment = True
        sqlpar = dict(lenMin = 5000)
        sql = """select name,count(*) as cnt from scaff_pred group by name order by cnt desc"""
        comment = "Count of assignments grouped by reference name"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog="\n")
        sql = """select name_species,count(*) as cnt from scaff_pred group by name_species order by cnt desc"""
        comment = "Count of assignments grouped by reference species name"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog="\n")
        sql = """select name_genus,count(*) as cnt from scaff_pred group by name_genus order by cnt desc"""
        comment = "Count of assignments grouped by reference genus name"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog="\n")
        sql = """select name_family,count(*) as cnt from scaff_pred group by name_family order by cnt desc"""
        comment = "Count of assignments grouped by reference family name"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment)
        outCsv.close()
        db.close()

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmClassifierApp)

