### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application for  building and using IMM-based taxonomic classifier in the spirit of Phymm"""

from MGT.ImmApp import *
from MGT.Taxa import *
from MGT.App import *
from MGT.DirStore import *
from MGT.SeqDbFasta import *
from MGT.ArchiveApp import *
from MGT.ImmClassifierBench import *

import functools

def fastaReaderFilterNucDegen(fastaReader,extraFilter=None):
    compr = SymbolRunsCompressor(sym="N",minLen=1)
    nonDegenSymb = "ATCGatcg"
    def line_gen():
        for rec in fastaReader.records():
            hdr = rec.header()
            seq = compr(rec.sequence())
            if not checkSaneAlphaHist(seq,nonDegenSymb,minNonDegenRatio=0.98):
                print "WARNING: ratio of degenerate symbols is too high, "+\
                        "skipping sequence with id %s" % (rec.getId(),)
            else:
                rec = extraFilter(hdr,seq)
                if rec:
                    hdr,seq = rec
                    yield hdr
                    for line in seqToLines(seq):
                        yield line
    return FastaReader(line_gen())

class TaxaPred(Struct):
    """Class to represent taxonomic predictions - one taxid per sample"""
    pass

class ImmClassifierApp(App):
    """App-derived class for building and using IMM-based taxonomic classifier in the spirit of Phymm.
    The main difference with Phymm is that we build IMMs for higher-level clades by pulling sequence
    data for subnodes.
    This class can be mostly viewed as imposing a TaxaTree structure onto ImmApp."""

    batchDepModes = ("predict","score","train","make-ref-seqdb",
            "fin-ref-seqdb","make-bench","bench-one-frag-len",
            "bench")

    ## Special taxonomy ID value to mark rejected samples 
    rejTaxid = 0
    
    # Flag set on a TaxaNode node to show that the model should not be
    # trained for this node
    TRAIN_SEL_STATUS_IGNORE = 0x00
    # Flag set on a TaxaNode node to show that the model can be trained
    # because the node has directly attached sequence
    TRAIN_SEL_STATUS_DIRECT = 0x01
    # Flag set on a TaxaNode node to show that the model can be trained
    # because the node has indirectly attached sequence
    TRAIN_SEL_STATUS_INDIRECT = 0x02

    appOptHelp = \
    """This is a unified driver program for ICM-based classification and model training.
    Quick start if you only want to classify with the pre-built default models:
    In Bash shell:

    source <MGT_HOME>/etc/mgtaxa.shrc
    python $MGT_HOME/bin/454_contig_read_cnt.py < 454ReadStatus.txt > weights.csv
    mgt-icm-classifier --mode predict --inp-seq 454LargeContigs.fna \\
            --inp-seq-attrib weights.csv --pred-min-len-samp 1000 \\
            --pred-out-dir my_results --run-mode batchDep \\
            --lrm-user-options "-P 0116"

    Where <MGT_HOME> should be replaced with MGTAXA installation directory. 
    In this example: 
    
    The local resource manager (LRM) is SGE and is passed user-specifc project code
    (replace it with your own and keep the proper quotation as above).

    All input sequences with length less than 1000 bp (--pred-min-len-samp) will
    be ignored.

    The '--run-mode batchDep' executes the program in a distributed mode (under a batch
    manager). You can switch to local sequential execution in a single process with
    --run-mode inproc.
    
    The weight file from 454 Newbler assembly output was generated with a provided helper
    script that extracts as weight the number of reads per contig. You can supply any other
    number as a per-sequence weight.
    That file, if provided, will be used when generating the aggregated clade abundance tables.
    If you are not working with assembly or not caring about the aggregated reports, skip
    the weight file generation and omit the --inp-seq-attrib option.
    
    It is required that the defline of your input FASTA file contained a unique ID for
    each sequence (the part between the '>' and the first blank character).

    The results will be in 'my_results' directory defined by --pred-out-dir options.

    In batchDep run-mode, the application currently leaves some service files in the starting
    directory after it finishes. You might want to create a temporary directory and start the
    program from where, so that it will be easy to clean up afterwards by removing that directory.

    Other examples:

    Build the sequence DB for model training from NCBI RefSeq multi-FASTA file(s).
    This currently filters the input by excluding plasmids and taxa w/o enough total sequence.
    
    mgt-icm-classifier --mode make-ref-seqdb \\
            --inp-ncbi-seq 'test_data/seqdb-fasta/*.fasta.gz' \\
            --db-seq tmp.db-seq --run-mode batchDep \\
            --lrm-user-options '-P 0413'
    
    Train models based on a sequence DB built by make-ref-seqdb step.

    mgt-icm-classifier --mode train --db-seq tmp.db-seq \\
            --db-imm tmp.imm --run-mode batchDep \\
            --lrm-user-options '-P 0413'

    Make a prediction for each sequence in the --inp-seq multi-FASTA file
    against the --db-imm database of models. The output results are stored in
    --pred-out-dir. Per-sequence predictions are stored in a CSV file. Aggregated
    counts per clade at various taxonomic levels are provided in the stats sub-directory,
    along with the auto-generated graphs. The same data is provided in a SQLite file.
    Several extra options can change default locations of the individual output files.
    If you omit the --db-imm option, the program will try to use a central DB of models
    configured for this installation.

    mgt-icm-classifier --mode predict --inp-seq 195.fasta.gz --db-imm tmp.imm \\
            --pred-min-len-samp 1000 --pred-out-dir tmp.results \\
            --run-mode batchDep --lrm-user-options '-P 0413'
    """

    @staticmethod
    def pathMultiOptToAbs(opt,name):
        val = opt[name]
        if val is not None:
            assert not isinstance(val,str),"Expected a list option here"
            opt[name] = [ os.path.abspath(x) for x in val ]

    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        
        optChoicesMode = ("make-ref-seqdb","train","predict","score","make-bench",
                "setup-train","proc-scores",
                "export-predictions","stats-pred",
                "extract-ref-seqdb","fin-ref-seqdb",
                "fin-ref-seqdb-batch",
                "make-bench-batch",
                "make-bench-fin",
                "proc-bench-scores",
                "bench-one-frag-len",
                "bench")
        
        option_list = [
            
            make_option("-m", "--mode",
            action="store", 
            type="choice",
            choices=optChoicesMode,
            dest="mode",
            help=(("What to do, choice of %s. The typical user entry points are: "+\
                "make-ref-seqdb,train,predict.") % (optChoicesMode,))),
            
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
                    "but only one entry - in training mode. If 'predict' "+\
                    "mode is used, the deafult will be the central database "+\
                    "under MGT_DATA."),
            
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
            dest="immIdToSeqIds",
            help="File that maps IMM IDs to lists of seq IDs during IMM training"),
            
            optParseMakeOption_Path(None, "--imm-ids",
            dest="immIds",
            help="File with list of IMM IDs to use in scoring and prediction. Default is all"+\
                    " IMMs from --imm-seq-ids"),
            
            optParseMakeOption_Path(None, "--inp-seq",
            dest="inpSeq",
            help="File with input FASTA sequence for prediction"),
            
            optParseMakeOption_Path(None, "--inp-seq-attrib",
            dest="sampAttrib",
            help="Optional tab-delimited file with extra attributes for each input sequence. "+\
                    "Currently must have two columns (no header row): sample id and weight. "+\
                    "Weight can be read count in a contig, and will be used when calculating "+\
                    "summary abundance tables."),

            optParseMakeOption_Path(None, "--inp-train-seq",
            dest="inpTrainSeq",
            help="File with input FASTA sequences for training extra user models"),
            
            optParseMakeOption_Path(None, "--inp-ncbi-seq",
            dest="inpNcbiSeq",
            help="File or shell glob with input NCBI FASTA sequences for training main set of models"),
            
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
            
            optParseMakeOption_Path(None, "--score-out-dir",
            dest="outDir",
            help="Directory name for output score files [--cwd/scores]"),
            
            optParseMakeOption_Path(None, "--out-score-comb",
            dest="outScoreComb",
            help="Output file for combined raw scores [--out-dir/combined.%s]"%(ImmApp.scoreSfx,)),
        
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
            dest="trainMinLenSamp",
            help="Min length of leaf node sequence to consider the node for training modesl: "+\
                    "default is 100000 if training custom models and 500000 if training "+\
                    "reference models. For training custom models, this means that sequences"+\
                    "shorter that this will be ignored"),
            
            make_option(None, "--train-max-len-samp-model",
            action="store", 
            type="int",
            default=10**9/2, #1G with rev-compl
            dest="trainMaxLenSampModel",
            help="Max length of sequence to use when training one model. If the available "+\
                    "sequence is longer, it will be subsampled [%default]"),
            
            make_option(None, "--new-tax-name-top",
            action="store", 
            type="string",
            default="mgt_reference",
            dest="newTaxNameTop",
            help="Root name for custom reference sequences"),
            
            optParseMakeOption_Path(None, "--pred-out-dir",
            dest="predOutDir",
            help="Output directory for classification results [--cwd/results]"),
            
            optParseMakeOption_Path(None, "--pred-out-taxa",
            dest="predOutTaxa",
            help="Output file with predicted taxa [--pred-out-dir/pred-taxa]"),
            
            optParseMakeOption_Path(None, "--pred-out-taxa-csv",
            dest="predOutTaxaCsv",
            help="Output CSV file with predicted taxa; default is --pred-out-taxa.csv"),
            
            optParseMakeOption_Path(None, "--trans-pred-out-taxa",
            dest="transPredOutTaxa",
            help="Existing output file with predicted taxa for custom training sequences to be used "+\
                    "in transitive classification [Optional]"),
            
            optParseMakeOption_Path(None, "--pred-out-stats-csv",
            dest="predOutStatsCsv",
            help="Output CSV file with statistics on predicted taxa [--pred-out-dir/stats/stats.csv]"),
            
            optParseMakeOption_Path(None, "--pred-out-stats-pdf",
            dest="predOutStatsPdf",
            help="Output PDF file with statistics on predicted taxa [--pred-out-taxa/stats/stats.pdf]"),
            
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
            
            make_option(None, "--incremental-work",
            action="store",
            type="int",
            default=0,
            dest="incrementalWork",
            help="Work incrementally wherever possible (restart mode). Currently, this will "+\
                    "affect 'train' and 'predict' modes [%default]. Set to non-zero number to "+\
                    "activate"),
            
            make_option(None, "--reduce-scores-early",
            action="store",
            type="int",
            default=1,
            dest="reduceScoresEarly",
            help="Try to perform score reduction/model selection as early as possible "+\
                    "and incrementally during the execution of the distributed scoring "+\
                    "pipeline in order to minimize per-process memory requirements and "+\
                    "overall speed. This is the default mode [set to 1]. It will be set to 0 "+\
                    "automatically during benchmark execution. If set to 0, a vector of all "+\
                    "model scores will be accumulated per each sample and processed in the "+\
                    "final reduction job-task."),
            
            make_option(None, "--bench-frag-len-list",
            action="store",
            type="string",
            default="100,400,1000,10000",
            dest="benchFragLenList",
            help="List of fragment lengths to generate and score in benchmark"),
            
            optParseMakeOption_Path(None, "--db-bench",
            dest="dbBench",
            help="BenchDb path, used either to create a new BenchDb from SeqDb "+\
                "or to use an existing BenchDb for benchmarking the models [-cwd/dbBench]"),
            
            optParseMakeOption_Path(None, "--db-bench-frag",
            dest="dbBenchFrag",
            help="Path of a single FASTA file that contains samples extracted from "+\
                "the BenchDb [--cwd/dbBenchFrag]"),
            
            make_option(None, "--db-bench-frag-len",
            action="store",
            type="int",
            default=400,
            dest="dbBenchFragLen",
            help="Length of fragments to generate in benchmark"),
            
            make_option(None, "--db-bench-frag-count-max",
            action="store",
            type="int",
            default=100,
            dest="dbBenchFragCountMax",
            help="Maximum count of fragments to select for benchmark per genome"),
            
            make_option(None, "--bench-n-lev-test-model-min",
            action="store",
            type="int",
            default=2,
            dest="benchNLevTestModelsMin",
            help="Minimum number of models below a node at a given taxonomic level required to evaluate "+\
                    "a benchmark performance for this node. The lowest possible number is two, because with "+\
                    "just one model we have no chance to produce the true prediction after we exclude one "+\
                    "model that contains the testing samples."),
            
            optParseMakeOption_Path(None, "--bench-out-dir",
            dest="benchOutDir",
            help="Output directory for benchmarking results [-cwd/benchResults]"),
            
            optParseMakeOption_Path(None, "--bench-out-csv",
            dest="benchOutCsv",
            help="Output CSV file with per-clade benchmarking results [--bench-out-dir/bench.csv]"),
            
            optParseMakeOption_Path(None, "--bench-out-aggr-csv",
            dest="benchOutAggrCsv",
            help="Output CSV file with aggregated benchmarking results [--bench-out-dir/bench.csv]"),
            
            optParseMakeOption_Path(None, "--bench-out-db-sqlite",
            dest="benchOutDbSqlite",
            help="Output SQLite database file with benchmarking results [--bench-out-dir/bench.sqlite]"),
        ]
        return Struct(usage = klass.appOptHelp+"\n"+\
                "Run with the --help options for a detailed description of individual arguments.",
                option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        globOpt = globals()["options"]
        if opt.isUndef("mode"):
            parser.error("--mode option is required")
        opt.setIfUndef("cwd",os.getcwd())
        if ( not opt.immDbArchive and not opt.immDb ):
            if opt.mode in ("predict","make-bench","bench"):
                opt.immDb = [ globOpt.icm.icmDb ]
        if isinstance(opt.benchFragLenList,str):
            opt.benchFragLenList = [ int(x) for x in opt.benchFragLenList.split(",") ]
        klass.pathMultiOptToAbs(opt,"immDb")           
        klass.pathMultiOptToAbs(opt,"immDbArchive")           
        opt.setIfUndef("immIdToSeqIds",pjoin(opt.cwd,"imm-seq-ids"))
        opt.setIfUndef("outDir",pjoin(opt.cwd,"scores"))
        opt.setIfUndef("predOutDir",pjoin(opt.cwd,"results"))
        opt.setIfUndef("seqDb",pjoin(opt.cwd,"seqDb"))
        opt.setIfUndef("dbBench",pjoin(opt.cwd,"dbBench"))
        opt.setIfUndef("dbBenchFrag",pjoin(opt.cwd,"dbBenchFrag.fna"))
        opt.setIfUndef("immIds",opt.immIdToSeqIds)
        opt.setIfUndef("outScoreComb",pjoin(opt.outDir,"combined"+ImmApp.scoreSfx))
        opt.setIfUndef("predOutTaxa",pjoin(opt.predOutDir,"pred-taxa"))
        opt.setIfUndef("predOutTaxaCsv",opt.predOutTaxa+".csv")
        opt.setIfUndef("predOutDbSqlite",opt.predOutTaxa+".sqlite")
        opt.setIfUndef("predOutStatsDir", pjoin(opt.predOutDir,"stats"))
        opt.setIfUndef("predOutStatsCsv", pjoin(opt.predOutStatsDir,"stats.csv"))
        opt.setIfUndef("predOutStatsPdf", pjoin(opt.predOutStatsDir,"stats.pdf"))
        opt.setIfUndef("newTaxidTop",mgtTaxidFirst)
        opt.setIfUndef("immDbWorkDir",pjoin(opt.cwd,"immDbWorkDir"))
        opt.setIfUndef("scoreWorkDir",pjoin(opt.cwd,"scoreWorkDir"))
        opt.setIfUndef("benchWorkDir",pjoin(opt.cwd,"benchWorkDir"))
        opt.setIfUndef("benchOutDir",pjoin(opt.cwd,"benchResults"))
        opt.setIfUndef("benchOutCsv",pjoin(opt.benchOutDir,"bench.csv"))
        opt.setIfUndef("benchOutAggrCsv",pjoin(opt.benchOutDir,"bench.aggr.csv"))
        opt.setIfUndef("benchOutDbSqlite",pjoin(opt.benchOutDir,"bench.sqlite"))
        if opt.predMode == "host": 
            opt.setIfUndef("rejectRanksHigher","order")
            opt.setIfUndef("predMinLenSamp",5000)
        elif opt.predMode == "taxa":
            opt.setIfUndef("rejectRanksHigher","superkingdom")
            opt.setIfUndef("predMinLenSamp",300)
        else:
            raise ValueError("Unknown --pred-mode value: %s" % (opt.predMode,))
        if opt.mode == "train":
            if opt.inpTrainSeq:
                defTrainMinLenSamp = 100000
            else:
                defTrainMinLenSamp = 500000
            opt.setIfUndef("trainMinLenSamp",defTrainMinLenSamp)
        if opt.mode == "make-ref-seqdb":
            globOpt = globals()["options"]
            opt.setIfUndef("inpNcbiSeq",pjoin(globOpt.refSeqDataDir,"microbial.genomic.fna.gz"))
        if opt.benchNLevTestModelsMin < 2:
            parser.error("--bench-n-lev-test-model-min must be higher than 1")

    
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
        makedir(opt.cwd)
        makedir(opt.immDbWorkDir)
        makedir(opt.scoreWorkDir)
        makedir(opt.benchWorkDir)
        #print "Options used:\n", opt
   

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
        elif opt.mode == "export-predictions":
            ret = self.exportPredictions(**kw)
        elif opt.mode == "proc-scores":
            ret = self.processImmScores(**kw)
        elif opt.mode == "setup-train":
            return self.setupTraining(**kw)
        elif opt.mode == "combine-scores":
            ret = self.combineScores(**kw)
        elif opt.mode == "make-ref-seqdb":
            ret = self.makeRefSeqDb(**kw)
        elif opt.mode == "extract-ref-seqdb":
            ret = self.extractRefSeqDb(**kw)
        elif opt.mode == "stats-pred":
            ret = self.statsPred(**kw)
        elif opt.mode == "fin-ref-seqdb":
            ret = self.finRefSeqDb(**kw)
        elif opt.mode == "fin-ref-seqdb-batch":
            ret = self.finRefSeqDbBatch(**kw)
        elif opt.mode == "make-bench":
            ret = self.makeBenchFromSeqDb(**kw)
        elif opt.mode == "make-bench-batch":
            ret = self.makeBenchFromSeqDbBatch(**kw)
        elif opt.mode == "make-bench-fin":
            ret = self.makeBenchFromSeqDbFin(**kw)
        elif opt.mode == "proc-bench-scores":
            ret = self.procBenchScores(**kw)
        elif opt.mode == "bench-one-frag-len":
            ret = self.benchOneFragLen(**kw)
        elif opt.mode == "bench":
            ret = self.bench(**kw)
        elif opt.mode == "bench-many-frag-len-fin":
            ret = self.benchManyFragLenFin(**kw)
        else:
            raise ValueError("Unknown mode value: %s" % (opt.mode,))
        return ret

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree(pklFile=self.opt.taxaTreePkl)
        return self.taxaTree

    def getTaxaLevels(self):
        if self.taxaLevels is None:
            #that assigns "level" and "idlevel" attributes to TaxaTree nodes,
            #when taxaTree is already loaded. Otherwise, you can use
            #taxaLevels.setTaxaTree() later.
            self.taxaLevels = TaxaLevels(self.taxaTree)
        return self.taxaLevels
    
    def getSeqDb(self):
        opt = self.opt
        if self.seqDb is None:
            self.seqDb = SeqDbFasta.open(opt.seqDb,mode="r") #"r"
            return self.seqDb

    def makeRefSeqDb(self,**kw):
        """Create reference SeqDb (the "main" SeqDb)"""
        opt = self.opt
        optI = copy(opt)
        optI.mode = "extract-ref-seqdb"
        app = self.factory(opt=optI)
        jobs = app.run(**kw)
        optI = copy(opt)
        optI.mode = "fin-ref-seqdb"
        app = self.factory(opt=optI)
        kwI = kw.copy()
        kwI["depend"] = jobs
        return app.run(**kwI)
    
    def extractRefSeqDb(self,**kw):
        """Extract reference SeqDb from input database sequences"""
        return self.makeNCBISeqDb(**kw)

    def makeNCBISeqDb(self,**kw):
        """Create SeqDb for training ICM model from NCBI RefSeq"""
        opt = self.opt
        self.seqDb = None
        filt = functools.partial(fastaReaderFilterNucDegen,
                extraFilter=lambda hdr,seq: None if "plasmid" in hdr.lower() else (hdr,seq) )
        seqDb = SeqDbFasta.open(path=opt.seqDb,mode="c")
        seqDb.setTaxaTree(self.getTaxaTree())
        seqDb.importByTaxa(glob.glob(opt.inpNcbiSeq),filt=filt)

        # For now, we filter by length when we build models.
        #taxids = seqDb.getTaxaList()
        #for taxid in taxids:
        #    taxidSeqLen = seqDb.seqLengths(taxid)["len"].sum()
        #    if taxidSeqLen < opt.trainMinLenSamp:
        #        seqDb.delById(taxid)

    def finRefSeqDb(self,**kw):
        """Finalize creation of SeqDb for training ICM model from NCBI RefSeq"""
        opt = self.opt
        seqDb = self.getSeqDb()
        taxids = seqDb.getIdList()
        taxids = n.asarray(taxids,dtype="O")
        nrnd.shuffle(taxids)
        jobs = []
        for taxidsBatch in n.array_split(taxids,min(1000,len(taxids))):
            optI = copy(opt)
            optI.mode = "fin-ref-seqdb-batch"
            optI.seqDbIds = taxidsBatch
            app = self.factory(opt=optI)
            jobs += app.run(**kw)
        return jobs
    
    def finRefSeqDbBatch(self,**kw):
        """Sub-task of finalizing creation of SeqDb for training ICM model from NCBI RefSeq"""
        opt = self.opt
        seqDb = self.getSeqDb()
        taxids = opt.seqDbIds
        for taxid in taxids:
            seqDb.finById(taxid)
            #print "DEBUG: finByid(%s) done" % (taxid,)

    ## Methods that generate benchmark dataset

    def makeBenchFromSeqDb(self,**kw):
        """Generate a benchmark dataset from SeqDb that was used to train ICM models"""
        opt = self.opt
        seqDb = self.getSeqDb()
        immDbs = [ ImmStore.open(path=immDb,mode='r') for immDb in opt.immDb ]
        bench = ImmClassifierBenchmark.open(path=opt.dbBench,mode="c")
        bench.seqDb = seqDb
        ids = bench.selectIdsDb(immDbs=immDbs)
        ids = n.asarray(ids,dtype="O")
        nrnd.shuffle(ids)
        jobs = []
        for idsBatch in n.array_split(ids,min(200,len(ids))):
            optI = copy(opt)
            optI.mode = "make-bench-batch"
            optI.benchIds = idsBatch
            app = self.factory(opt=optI)
            jobs += app.run(**kw)
        
        coOpt = copy(opt)
        coOpt.mode = "make-bench-fin"
        coOpt.benchIds = ids
        coApp = self.factory(opt=coOpt)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = coApp.run(**kw)
        return jobs
    
    def makeBenchFromSeqDbBatch(self,**kw):
        """Sub-task of benchmark generation"""
        opt = self.opt
        seqDb = self.getSeqDb()
        ids = opt.benchIds
        bench = ImmClassifierBenchmark.open(path=opt.dbBench,mode="r")
        #bench.save() should not be called with seqDb attribute assigned
        #or it will try to dump it as well
        bench.seqDb = seqDb
        for id in ids:
            bench.makeSample(idDb=id,
                    fragLen=opt.dbBenchFragLen,
                    fragCountMax=opt.dbBenchFragCountMax)
        
    def makeBenchFromSeqDbFin(self,**kw):
        """Final task of benchmark generation"""
        opt = self.opt
        ids = opt.benchIds
        bench = ImmClassifierBenchmark.open(path=opt.dbBench,mode="r")
        bench.catSamples(outFile=opt.dbBenchFrag,idsDb=ids)

    ## Methods that integrate all steps of benchmark generation and evaluation

    def bench(self,**kw):
        """Generate and evaluate the entire benchmark suite"""
        return self.benchManyFragLen(**kw)

    def benchManyFragLen(self,**kw):
        """Generate and evaluate benchmark for a range of fragment lengths."""
        opt = self.opt
        jobs = []
        optsI = []
        for fragLen in opt.benchFragLenList:
            optI = Struct()
            optI.runMode = opt.runMode
            optI.lrmUserOptions = opt.lrmUserOptions
            optI.mode = "bench-one-frag-len"
            optI.immDb = opt.immDb
            optI.seqDb = opt.seqDb
            optI.dbBenchFragLen = fragLen
            optI.dbBenchFragCountMax = opt.dbBenchFragCountMax
            optI.cwd = pjoin(opt.benchWorkDir,str(fragLen))
            #this sets default and derived options
            #@todo we need to create two kinds of properties in Struct():
            #fixed values and lazy-evaluated values (similar to variables
            #in "make" assigned with "=" and ":="). This will allow changing
            #just one option in a copy of the parent option object (e.g.
            #the working directory), with dependent options (e.g. various
            #filenames) dynamically recomputed.
            app = self.factory(opt=optI)
            optsI.append(app.opt)
            jobs += app.run(**kw)
        
        optI = copy(opt)
        optI.mode = "bench-many-frag-len-fin"
        optI.optsOneFragLen = optsI
        app = self.factory(opt=optI)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = app.run(**kw)

        return jobs
    
    def benchManyFragLenFin(self,**kw):
        """Combiner phase of benchManyFragLen()"""
        return kw.get("depend",None)
    
    def benchOneFragLen(self,**kw):
        """Generate and evaluate benchmark for a given fragment length."""
        opt = self.opt
        optI = copy(opt)
        optI.mode = "make-bench"
        app = self.factory(opt=optI)
        jobs = app.run(**kw)
        
        optI = copy(opt)
        optI.mode = "predict"
        optI.inpSeq = opt.dbBenchFrag
        optI.predMinLenSamp = 1
        optI.reduceScoresEarly = 0
        app = self.factory(opt=optI)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = app.run(**kw)
        
        optI = copy(opt)
        optI.mode = "proc-bench-scores"
        app = self.factory(opt=optI)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = app.run(**kw)

        return jobs

    ## Methods that assign training sequences to higher-level nodes

    def mapSeqToTree(self):
        """Assign list of SeqDB IDs to corresponding taxonomy tree nodes.
        In the current SeqDB version, each ID is a unique taxonomy id, so
        a one-element list will be assigned.
        @post Attribute 'leafSeqDbIds' is assigned to EVERY node and contains a list of IDs, possibly empty.
        The empty list will be a reference to a single shared object, so it should be treated as immutable"""
        opt = self.opt
        taxaTree = self.getTaxaTree()
        seqDb = self.getSeqDb()
        taxaList = [ taxid for taxid in seqDb.getTaxaList() \
                if seqDb.seqLengths(taxid)["len"].sum() >= opt.trainMinLenSamp ]
        emptyList = []
        taxaTree.setAttribute("leafSeqDbIds",emptyList,doCopy=False)
        for taxid in taxaList:
            taxaTree.getNode(taxid).leafSeqDbIds = [ taxid ]


    def pickSeqOnTree(self,maxSeqIdCnt):
        """Assign to each node of the taxonomy tree the list of SeqDB IDs in the subtree under that node.
        It picks a median count of IDs from each child node and itself, clipped by maxSeqIdCnt.
        """
        taxaTree = self.getTaxaTree()

        def actor(node):
            if node.isLeaf() and len(node.leafSeqDbIds)>0:
                node.pickedSeqDbIds = node.leafSeqDbIds
                node.trainSelStatus = self.TRAIN_SEL_STATUS_DIRECT
            else:
                chSeqIds = [c.pickedSeqDbIds for c in node.getChildren() if hasattr(c,"pickedSeqDbIds")]
                chSeqIds.append(node.leafSeqDbIds)
                chSeqIds = [ x for x in chSeqIds if len(x) > 0 ]
                if len(chSeqIds) > 0:
                    targLen = int(n.median([ len(x) for x in chSeqIds ]))
                    targLen = min(maxSeqIdCnt,targLen)
                    chSeqIdsSel = [ sampleBoundWOR(x,targLen) for x in chSeqIds ]
                    node.pickedSeqDbIds = sum(chSeqIdsSel,[])
                    assert len(node.pickedSeqDbIds) > 0
                    #We only train for nodes with directly attached sequence,
                    #or at least two subnodes with some (possibly indirectly
                    #attached) sequence. This is done to avoid creating a lineage
                    #of identical models when there is only one sequence sample
                    #at the bottom of the lineage. This implicitely selects the
                    #most specific assignment for such lineages, unless we change
                    #it at the prediction stage by looking at the location of models
                    #on the tree.
                    if len(node.leafSeqDbIds)>0:
                        node.trainSelStatus = self.TRAIN_SEL_STATUS_DIRECT
                    elif len(chSeqIds) >= 2:
                        node.trainSelStatus = self.TRAIN_SEL_STATUS_INDIRECT


        taxaTree.visitDepthBottom(actor)

    def defineImms(self):
        taxaTree = self.getTaxaTree()
        micNodes = [ taxaTree.getNode(taxid) for taxid in micTaxids ]
        cellNode = taxaTree.getNode(cellTaxid)
        #DEBUG:
        #cntLeafSeq = sum([ len(node.leafSeqDbIds)>0 for node in taxaTree.iterDepthTop() if node.isUnder(cellNode) ])
        #cntDirect = sum([ hasattr(node,"trainSelStatus") and node.trainSelStatus == self.TRAIN_SEL_STATUS_DIRECT for node in taxaTree.iterDepthTop() if node.isUnder(cellNode) ])
        immIdToSeqIds = {}
        ##@todo Make it controlled by a node selection algebra passed by the user as json expression
        for node in taxaTree.iterDepthTop():
            if hasattr(node,"trainSelStatus") and node.trainSelStatus != self.TRAIN_SEL_STATUS_IGNORE:
                doPick = False
                if node.isUnderAny(micNodes):
                    doPick = True
                elif node.isUnder(cellNode):
                    if node.trainSelStatus == self.TRAIN_SEL_STATUS_DIRECT:
                        doPick = True
                if doPick:
                    #DEBUG:
                    #print "Training for: ", node.lineageStr()
                    immIdToSeqIds[node.id] = node.pickedSeqDbIds
        dumpObj(immIdToSeqIds,self.opt.immIdToSeqIds)

    def _customTrainSeqIdToTaxaName(self,seqid):
        """Generate a name for new TaxaTree node from sequence ID.
        To be used both when defining tree nodes for custom training sequences,
        as well as when joining with predictions for these sequences for
        transitive annotation"""
        return self.opt.newTaxNameTop+"_"+seqid

    def _customTrainSeqIdToTaxid(self,seqid):
        """Generate a name for new TaxaTree node from sequence ID.
        To be used both when joining with predictions for these sequences for
        transitive annotation.
        @return TaxaNode ID or None if seqid not found"""
        taxaName = self._customTrainSeqIdToTaxaName(seqid)
        node = self.getTaxaTree().searchName(taxaName) #returns a list
        if len(node) > 1:
            raise ValueError("Custom taxa name must be unique in the tree, found: %s" % \
                ",    ".join(["%s" % _n for _n in node]))
        elif len(node) == 1:
            return node[0].id
        else:
            return None

    def _customTrainTaxidToSeqId(self,taxid):
        """Return seq id corresponding to a custom TaxaTree node.
        To be used both when joining with predictions for these sequences for
        transitive annotation.
        The implementation must match _customTrainSeqIdToTaxaName()"""
        taxaTree = self.getTaxaTree()
        node = taxaTree.getNode(taxid)
        try:
            return node.name.split(self.opt.newTaxNameTop+"_")[1]
        except IndexError:
            return None

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
                taxaNode = TaxaNode(id=nextNewTaxid,name=self._customTrainSeqIdToTaxaName(seqid),rank=unclassRank,divid=dividEnv,names=list())
                taxaNode.setParent(newTaxaTop)
                # that will be used by the following call to defineImms()
                taxaNode.pickedSeqDbIds = [ taxaNode.id ]
                taxaNode.trainSelStatus = self.TRAIN_SEL_STATUS_DIRECT
                fastaWriter = seqDb.fastaWriter(id=nextNewTaxid,lineLen=80)
                fastaWriter.record(header=hdr,sequence=seq)
                fastaWriter.close()
                seqDb.finById(id=nextNewTaxid)
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
        makeFilePath(opt.outScoreComb)
        immScores = openImmScores(opt,fileName=opt.outScoreComb,mode="w")
        immScores.catImms(fileNames=opt.outSubScores)
        immScores.close()

    def predict(self,**kw):
        """Score input sequences and predict taxonomy"""
        opt = self.opt
        
        optI = copy(opt)
        optI.mode = "score"
        app = self.factory(opt=optI)
        jobs = app.run(**kw)
        
        optI = copy(opt)
        optI.mode = "proc-scores"
        app = self.factory(opt=optI)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = app.run(**kw)
        return jobs

    def _maskScoresNonSubtrees(self,taxaTree,immScores,posRoots):
        """Set to a negative infinity (numpy.NINF) all columns in score matrix that point to NOT subtrees of posRoots nodes.
        This is used to mask all scores pointing to other than bacteria, archaea or eukaryots."""
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

    def _roundUpPredictions(self,predTaxids,rootNodes,rank):
        """Make predictions less specific for selected subtrees.
        For example, you can use this to make each eukaryotic prediction no
        more specific than a phylum (or a nearest defined rank).
        @param taxaTree TaxaTree instance
        @param predTaxids an array of predicted taxids, one per sample
        @param roots a sequence of TaxaNode instances defining subtrees that
        will be affected
        @param rank lowest rank allowed for selected subtrees"""
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        rankPos = taxaLevels.getLevelPos()[rank]
        for (i,taxid) in enumerate(predTaxids):
            if not taxid == self.rejTaxid:
                node = taxaTree.getNode(taxid)
                if node.isSubnodeAny(rootNodes):
                    lin = taxaLevels.lineageFixedList(node,null=None,format="node",fill="up-down")
                    rankNode = lin[rankPos]
                    #assert rankNode is not None,("Logic error for predNode=%s" % (predNode,))
                    #condition below should probably always be true...
                    if rankNode:
                        predTaxids[i] = rankNode.id

    def scoresToPredictions(self,immScores,**kw):
        """Generate taxonomic predictions from raw IMM scores.
        This handles both classification of bacterial sequences and host assignment 
        for viral sequences.
        @param immScores ImmScores object.
        
        Other parameters are taken from self.opt.
        @param rndScoreComb File with ImmScores object for random query sequences
        @param predOutTaxa Output file with predicted taxa
        """
        sc = immScores.getData()
        #assume idImm are str(taxids):
        idImm = n.asarray(sc.idImm[:],dtype=int)
        kind = immScores.getKind()
        if kind == "ImmScoresReduced":
            predTaxids = idImm
        elif kind == "ImmScoresDenseMatrix":
            #scRnd = loadObj(opt.rndScoreComb)
            #sc.score = self._normalizeScores(sc,scRnd)
            #normalize to Z-score along each row
            #sc.score = ((sc.score.T - sc.score.mean(1))/sc.score.std(1)).T
            #normalize to Z-score over entire matrix
            #sc.score = ((sc.score - sc.score.mean())/sc.score.std())
            score = sc.score[:]
            #taxaTree = self.getTaxaTree()
            #taxaLevels = self.getTaxaLevels()
            #micRoots = [ taxaTree.getNode(taxid) for taxid in micTaxids ]
            #virRoot = taxaTree.getNode(virTaxid)
            #cellRoot = taxaTree.getNode(cellTaxid)
            #scVirRoot = score[:,idImms == virTaxid][:,0]
            #scCellRoot = score[:,idImms == cellTaxid][:,0]
            #cellTopColInd = n.concatenate([ n.where(idImms == taxid)[0] for taxid in micTaxids ])
            #scCellRootMax = score[:,cellTopColInd].max(1)
            #self._maskScoresNonSubtrees(taxaTree,immScores=sc,posRoots=(cellRoot,))
            #if opt.rejectRanksHigher is not "superkingdom":
            #    self._maskScoresByRanks(taxaTree,immScores=sc)
            #topTaxids = self._taxaTopScores(taxaTree,immScores=sc,topScoreN=10)
            argmaxSc = score.argmax(0)
            #maxSc = score.max(0)
            predTaxids = idImm[argmaxSc]
        else:
            raise ValueError("Unknown kind of score object: %s" % (kind,))
        return predTaxids

    def filterPredictions(self,immScores,predTaxids,**kw):
        """Generate taxonomic predictions from raw IMM scores.
        This handles both classification of bacterial sequences and host assignment 
        for viral sequences.
        @param immScores ImmScores object.
        @param predTaxids[in,out] array with predicted taxonomy per sample
        
        Other parameters are taken from self.opt.
        """
        opt = self.opt
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        micRoots = [ taxaTree.getNode(taxid) for taxid in micTaxids ]
        virRoot = taxaTree.getNode(virTaxid)
        cellRoot = taxaTree.getNode(cellTaxid)
        #predTaxids[maxSc<0] = self.rejTaxid
        # this will reject any sample that has top level viral score more
        # that top level cellular org score, on the assumption that easily
        # assignbale viruses should be closer to cellular orgs than to other
        # viruses.
        # Result: removed 90 out of 250 samples, with no change in specificity.
        #predTaxids[scVirRoot>scCellRoot] = self.rejTaxid
        #predTaxids[scVirRoot>scCellRootMax] = self.rejTaxid
        # this excluded 10 out of 430, no change in specificity
        #predTaxids[scVirRoot>=maxSc] = self.rejTaxid
        
        # This is not the same as _maskScoresByRanks() above (because this
        # will actually assign a reject label).
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

        #Round-up euk predictions to the phylum level because of high sequence identities between
        #lower order clades
        self._roundUpPredictions(predTaxids,rootNodes=(taxaTree.getNode(eukTaxid),),rank="phylum")
    
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
        immScores = openImmScores(opt,fileName=opt.outScoreComb,mode="r")
        predTaxids = self.scoresToPredictions(immScores=immScores,**kw)
        self.filterPredictions(immScores=immScores,predTaxids=predTaxids)
        sc = immScores.getData()
        pred = TaxaPred(idSamp=sc.idSamp[:],predTaxid=predTaxids,lenSamp=sc.lenSamp[:])
        makeFilePath(opt.predOutTaxa)
        dumpObj(pred,opt.predOutTaxa)
        self.exportPredictions()
        #todo: "root" tax level; matrix of predhost idlevel's; batch imm.score() jobs; score hist;

    def exportPredictions(self,**kw):
        """
        Export predictions and summary statistics.
        Parameters are taken from self.opt.
        @param predOutTaxa Output file with predicted taxa
        @param sampAttrib Optional input tab-delimited file with extra per-sample
        attributes. Currently it should have two columns: sample id and 
        weight. Weight can be a number of reads in a given contig if
        samples are assembly contigs. The clade counts will be multiplied
        by these weights when aggregate summary tables are generated.
        """
        opt = self.opt
        self._exportPredictions()
        self._reExportPredictionsWithSql()
        self.statsPred()

    def _buildCustomTaxidToRefTaxidMap(self,predOutTaxa):
        """Helper method for transitive classification.
        This loads a prediction file created in a separate
        run for assigning custom model sequences to a reference DB,
        and returns a dict that maps custom taxonomic id generated
        here to assigned reference DB taxid"""
        opt = self.opt
        taxaPred = loadObj(predOutTaxa)
        custToRef = {}
        for idS,taxidRef,lenS in it.izip(taxaPred.idSamp,taxaPred.predTaxid,taxaPred.lenSamp):
            taxidCust = self._customTrainSeqIdToTaxid(idS)
            if taxidCust:
                custToRef[taxidCust] = taxidRef
        return custToRef

    def _exportPredictions(self):
        opt = self.opt    
        taxaPred = loadObj(opt.predOutTaxa)
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        levNames = taxaLevels.getLevelNames("ascend")
        if opt.transPredOutTaxa:
            transTaxaMap = self._buildCustomTaxidToRefTaxidMap(opt.transPredOutTaxa)
        else:
            transTaxaMap = None
        makeFilePath(opt.predOutTaxaCsv)
        out = openCompressed(opt.predOutTaxaCsv,"w")
        flds = ["id","len","taxid","name","rank"]
        for levName in levNames:
            flds += ["taxid_"+levName,"name_"+levName]
        w = csv.DictWriter(out, fieldnames=flds, restval='null',dialect='excel-tab')
        w.writerow(dict([(fld,fld) for fld in flds]))
        predPerSample = 1
        if predPerSample > 1:
            predTaxidRows = taxaPred.topTaxid[:,:predPerSample]
        else:
            predTaxidRows = taxaPred.predTaxid[:,n.newaxis]
        for idS,taxids,lenS in it.izip(taxaPred.idSamp,predTaxidRows,taxaPred.lenSamp):
            if lenS < opt.predMinLenSamp:
                continue
            for taxid in taxids:
                row = dict(id=idS,len=lenS)
                if transTaxaMap:
                    if taxid != self.rejTaxid:
                        taxid = transTaxaMap[taxid]
                if taxid != self.rejTaxid:
                    node = taxaTree.getNode(taxid)
                    row["name"] = node.name
                    row["rank"] = node.linn_level
                    lin = taxaLevels.lineageFixedList(node,null=None,format="node",fill="up-down")
                    for (levName,levNode) in it.izip(levNames,lin):
                        if levNode:
                            row["name_"+levName] = levNode.name
                            row["taxid_"+levName] = levNode.id
                row["taxid"] = taxid
                w.writerow(row)
        out.close()
        
    def _reExportPredictionsWithSql(self):
        opt = self.opt    
        makeFilePath(opt.predOutDbSqlite)
        db = DbSqlLite(dbpath=opt.predOutDbSqlite,strategy="exclusive_unsafe")
        db.createTableFromCsv(name="scaff_pred_1",
                csvFile=opt.predOutTaxaCsv,
                hasHeader=True,
                fieldsMap={"len":SqlField(type="integer")},
                indices={"names":("id",)})

        if opt.sampAttrib:
            db.createTableFromCsv(name="scaff_attr",
                    csvFile=opt.sampAttrib,
                    hasHeader=False,
                    fieldsMap={0:SqlField(name="id"),
                        1:SqlField(name="weight",type="real")},
                    indices={"names":("id",)})
        else:
            db.createTableAs("scaff_attr",
                    """\
                    select distinct id,1.0 as weight
                    from scaff_pred_1""",
                    indices={"names":("id",)})
        #make sure all samples can be matched up with attribute records
        db.executeAndAssertZero(\
                """select count(*) from scaff_pred_1 a
                where a.id not in (select id from scaff_attr)
                limit 1
                """,
                message="Some samples do not have matching sample attribute records - "+\
                        "check your --inp-seq-attrib file - it has to be strictly tab-delimited "+\
                        "with no extra spaces around fields.")
        sql = \
                """select a.*,b.weight as weight
                from scaff_pred_1 a, scaff_attr b
                where a.id = b.id
                """
        db.createTableAs("scaff_pred",sql,indices={"names":("id","taxid","name","len")})
        db.exportAsCsv("select * from scaff_pred",
            opt.predOutTaxaCsv,
            comment=None,
            sqlAsComment=False)
        db.close()

    def _sqlReport(self,db,dbTable,levName,outCsv=None):
        if levName:
            fldGrp = "name_"+levName
        else:
            fldGrp = "name"
        dbTableRep = dbTable+"_grp_"+levName
        sql = """\
                select %(fldGrp)s as clade,
                sum(weight) as sum_weight,
                count(*) as cnt_samp,
                sum(len) as len_samp,
                avg(len) as avg_len_samp 
                from %(dbTable)s
                group by %(fldGrp)s
                order by sum_weight desc
        """ % dict(fldGrp=fldGrp,dbTable=dbTable)
        db.createTableAs(dbTableRep,sql)
        if outCsv:
            comment = "Count of assignments grouped by %s" % \
                    (levName if levName else "lowest assigned clade",)
            db.exportAsCsv("""\
                    select * from %(dbTableRep)s
                    order by sum_weight desc
                    """ % dict(dbTableRep=dbTableRep),
                    outCsv,
                    comment=comment,
                    sqlAsComment=False,
                    epilog="\n")
        return dbTableRep
    
    def _graphicsReport(self,db,dbTableRep,levName,outBackend,fldRep="sum_weight",maxClades=20):
        import matplotlib
        matplotlib.use('AGG')
        from MGT import Graphics
        data=db.selectAll("""select
            clade, %(fldRep)s            
            from %(dbTableRep)s 
            order by %(fldRep)s desc
            limit %(maxClades)s
            """ % (dict(dbTableRep=dbTableRep,maxClades=maxClades,fldRep=fldRep)))
        Graphics.barHorizArea(data=data,
                xLabel="Count of assignments",
                yLabel=("Assigned %s" % (levName,)) if levName else "Lowest assigned clade",
                outBackend=outBackend)


    def statsPred(self,**kw):
        """Create aggregate tables,figures and csv files to show statistics of predictions"""
        opt = self.opt
        taxaLevels = self.getTaxaLevels()
        levNames = taxaLevels.getLevelNames("ascend")
        db = DbSqlLite(dbpath=opt.predOutDbSqlite,strategy="exclusive_unsafe")
        rmrf(opt.predOutStatsDir)
        makeFilePath(opt.predOutStatsCsv)
        outCsv = openCompressed(opt.predOutStatsCsv,'w')
        from matplotlib.backends.backend_pdf import PdfPages
        outBackend = PdfPages(opt.predOutStatsPdf)
        sqlAsComment = True
        sqlpar = dict(lenMin = opt.predMinLenSamp)
        db.ddl("""\
                create temporary view scaff_pred_filt as 
                select * from scaff_pred 
                where len >= %(lenMin)s""" % sqlpar)
        for levName in levNames+[""]: #"" is for lowest clade
            dbTableRep = self._sqlReport(db=db,dbTable="scaff_pred_filt",levName=levName,outCsv=outCsv)
            self._graphicsReport(db=db,dbTableRep=dbTableRep,levName=levName,fldRep="sum_weight",
                    outBackend=outBackend,maxClades=20)
        outCsv.close()
        outBackend.close()
        db.close()

    def procBenchScores(self,**kw):
        """Process benchmark IMM scores and generate performance metrics.
        Parameters are taken from self.opt.
        @param outScoreComb File with ImmScores object
        @param benchOutCsv Output file with per clade metrics in CSV format
        @param benchOutAggrCsv Output file with aggregated metrics in CSV format
        @param benchOutDbSqlite Output file with metrics in SQLite format
        """
        opt = self.opt
        makeFilePath(opt.benchOutCsv)
        makeFilePath(opt.benchOutDbSqlite)
        db = DbSqlLite(dbpath=opt.benchOutDbSqlite,strategy="exclusive_unsafe")
        storesDb = self.prepBenchSql(db=db)
        self.procBenchScoreMatr(db=db)
        sqlBenchMetr = ImmClassifierBenchMetricsSql(db=db,
                taxaLevelsTbl=storesDb.storeDbLev.tblLevels,
                taxaNamesTbl=storesDb.storeDbNode.tblNames)
        sqlBenchMetr.makeMetrics(nLevTestModelsMin=opt.benchNLevTestModelsMin,
                csvAggrOut=opt.benchOutAggrCsv)

    def prepBenchSql(self,db):
        """Prepare lookup tables in benchmark SQL database.
        @param db Sql object"""
        opt = self.opt
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        print "DEBUG: Loaded taxa tree and level"
        tblSfx = ""
        storeDbNode = NodeStorageDb(db=db,tableSfx=tableSfx)
        storeDbNode.saveName(taxaTree)
        storeDbLev = LevelsStorageDb(db=db,tableSfx=tableSfx)
        storeDbLev.save(taxaLevels)
        return Struct(storeDbLev=storeDbLev,storeDbNode=storeDbNode)

    def procBenchScoreMatr(self,db,**kw):
        """Process benchmark IMM scores and generate (test,pred) pairs in SQL.
        @param db Sql object

        Other parameters are taken from self.opt.
        @param outScoreComb File with ImmScores object
        @param benchOutCsv Output file with metrics in CSV format
        @param benchOutDbSqlite Output file with metrics in SQLite format
        """
        opt = self.opt
        immScores = openImmScores(opt,fileName=opt.outScoreComb,mode="r")
        sc = immScores.getData()
        #assume idImm are str(taxids):
        idImm = n.asarray(sc.idImm[:],dtype=int)
        kind = immScores.getKind()
        assert kind == "ImmScoresDenseMatrix","Matrix Score dataset is required for benchmark"
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        print "DEBUG: Loaded taxa tree and level"
        eukRoot = taxaTree.getNode(eukTaxid)
        levNames = taxaLevels.getLevelNames("ascend")
        levIdNames = taxaLevels.getLevelIdForNames(levNames)
        levPosSuperking = levNames.index("superkingdom")
        #this value as prediction result should be treated as
        #always incorrect when computing performance metrics
        wrongTaxid = self.rejTaxid - 1

        def _mapModelIds(taxaTree,idImm):
            """Map model ids to the taxonomy tree and store for each node the
            number of immediate children that have models in their subtrees.
            This will be used to include into the performance counters only
            those samples for which a correct model exists at a given exclusion
            level"""
            root = taxaTree.getRootNode()
            root.setAttribute("tmpHasMod",0)
            for taxid in idImm:
                node = taxaTree.getNode(taxid)
                node.tmpHasMod = 1
            root.setTotal("tmpHasMod","tmpCntMod")
            root.setReduction(extractor=lambda n: 0,
                    dstAttr="tmpChWMod",
                    childExtractor=lambda n: n.tmpCntMod > 0)


        _mapModelIds(taxaTree=taxaTree,idImm=idImm)
        print "DEBUG: mapped models to tree"


        create table  
        ( 
        db.ddl("""
        create table bench_samp 
        ( 
        i_samp integer NOT NULL, 
        i_lev_exc integer NOT NULL,
        i_lev_per integer NOT NULL,
        i_lev_test_real integer,
        i_lev_pred_real integer,
        taxid_lev_test integer,
        taxid_lev_pred integer,
        taxid_superking integer,
        n_lev_test_models integer
        )
        """,dropList=["table bench_samp"])
        inserter = db.makeBulkInserter(table="bench_samp")
        def _outTestRec(iS,iLevFS,iLevFP,nodLinFS_FP,nodLinFP,nodSuperkingS,
                inserter=inserter,levIdNames=levIdNames,
                missTaxid=missTaxid):
            rec = (iS,
                    levIdNames[iLevFS],
                    levIdNames[iLevFP],
                    nodLinFS_FP.idlevel,
                    nodLinFP.idlevel if nodLinFP is not None else levIdNames[iLevFP],
                    nodLinFS_FP.id,
                    nodLinFP.id if nodLinFP is not None else wrongTaxid,
                    nodSuperkingS.id,
                    nodLinFS_FP.tmpChWMod)
            inserter(rec)
        sampTaxids = set()
        for iS,idS in enumerate(sc.idSamp[:]):
            taxid,itS = idS.strip().split('_')
            taxid,itS = int(taxid),int(itS)
            sampTaxids.add(taxid)
        #If we select columns directly from hdf file with row-optimized
        #chunkshape, the performance is horrible and we spend 100% CPU in
        #system mode calls. For now, we load entire array into numpy.
        #@todo add copy to column-optimized chunkshape to make the final 
        #hdf5 (per pytables' author recipe), then select from file directly.
        score = sc.score[:]
        print "DEBUG: Loaded score matrix of total elements %i" % (score.shape[0]*score.shape[1],)
        #rows are models, columns are samples
        #samples ID was generate by bench generator and is "$taxid_$sampIndInTaxid"
        #we cache variables that depend only on true sample taxid or predicted sample
        #taxid because the former go in series, and the later often too (when predictions
        #are correct, for instance)
        tidS = None #taxid sample (true value)
        nodS = None #node(taxid sample)
        linS = None #lineage sample
        linSetS = None #lineage sample as set
        linFS = None #lineage fixed list for sample
        tidP = None #taxid predicted
        nodP = None #node(taxid predicted)
        linP = None #lineage predicted
        fill = None #None #"up-down" "up"
        assert fill is None,"The consequences of using fill for fixed lineage are not fully understood yet"
        iSel = sorted(sampleRatioWOR(n.arange(len(sc.idSamp)).tolist(),0.1))
        for iS,idS in it.izip(iSel,sc.idSamp[iSel]):
        #for iS,idS in enumerate(sc.idSamp):
            taxid,itS = idS.strip().split('_')
            taxid,itS = int(taxid),int(itS)
            if tidS is None or tidS != taxid:
                tidS = taxid
                nodS = taxaTree.getNode(tidS)
                linS = nodS.lineage()
                linSetS = set(linS)
                linFS = taxaLevels.lineageFixedList(nodS,null=None,format="node",fill=fill)
            ordP = (-(score[:,iS])).argsort()
            startOrdP = 0
            nodSuperkingS = linFS[levPosSuperking]
            #nodLinFS_Prev = None
            #iLevFS is exclusion level
            for iLevFS,nodLinFS in enumerate(linFS[:-1]):
                if nodLinFS is not None: 
                    #The condition above checks that the exclusion level is present in lineage 
                    #of a true node. Clearly, we cannot speak of excluding related models 
                    #at a certain level if the test sample does not have that level defined.
                    #So, the above omits such records from consideration regardless of the
                    #predicted taxa.
                    for iOrdP,vOrdP in enumerate(ordP,startOrdP):
                        tidP = idImm[vOrdP]
                        nodP = taxaTree.getNode(tidP)
                        #first condition excludes mixed models
                        #TMP:tidP in sampTaxids and
                        #not nodP.isUnder(eukRoot) and
                        if \
                            ((nodP not in linSetS \
                            and not nodP.isUnder(nodLinFS))):
                                linFP = taxaLevels.lineageFixedList(nodP,null=None,format="node",fill=fill)
                                #climb the levels and output (test,pred) pairs
                                #iLevFP is prediction level > iLevFS
                                for iLevFP,nodLinFP in enumerate(linFP[iLevFS+1:],iLevFS+1):
                                    nodLinFS_FP = linFS[iLevFP] #true node at prediction level
                                    if nodLinFS_FP is not None: 
                                        #The above checks that the prediction level is 
                                        #present in lineage of a true (test) node
                                        #Assuming the "fill" parameter for the fixed
                                        #lineage generation is None (so no fill), if the predicted
                                        #node does not have that level in its lineage (nodLinFP is None),
                                        #then it is a different lineage at that level and, therefore,
                                        #an incorrect prediction. We save a special value of "always
                                        #wrong taxid" as prediction value into a (test,predict) 
                                        #database record.
                                        #This reasoning is not that obvious when the fill is used
                                        #(although should probably work when fill=="up"), so we
                                        #restrict fill to None for now.

                                        #Check that sample is not the only child with model at this 
                                        #prediction level.
                                        #Phymm Methods for unexplained reason says they required 
                                        #"two others", which must mean 3 total, so we also save the 
                                        #number as part of the record to filter on it later when 
                                        #computing the metrics.
                                        if nodLinFS_FP.tmpChWMod >= 2: 
                                            _outTestRec(iS,iLevFS,iLevFP,nodLinFS_FP,nodLinFP,nodSuperkingS)
                                break
                    #next exclusion level is above current, its exclusion zone contains 
                    #the current zone,
                    #so it can ignore already scanned and rejected model indices
                    #and start from the last accepted index (careful that it does 
                    #not depend on the presence of levels in the lineage, but only
                    #on sub-tree relationship (it does not now).
                    #startOrdP = iOrdP
                    #nodLinFS_Prev = nodLinFS
            del ordP
        inserter.close()
        del score
        immScores.close()
        db.createIndices(table="bench_samp",
        names=[
        "i_lev_exc",
        "i_lev_per",
        "taxid_lev_test",
        "taxid_lev_pred",
        "taxid_superking",
        "n_lev_test_models"
        ])

        db.close()

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmClassifierApp)

