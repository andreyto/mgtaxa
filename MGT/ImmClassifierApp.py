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
from MGT.ArchiveGpgApp import *
from MGT.ArchiveApp import *
from MGT.ImmClassifierBench import *
from MGT.ImmClassifierAppUtils import *
from MGT.GraphicsKrona import *
from MGT.ImmScoreRescaler import *
from MGT.TaxaPred import TaxaPred
from MGT.FastaSplitters import *
from MGT.SeqDbFastaApp import *

import functools



class ImmClassifierApp(App):
    """App-derived class for building and using IMM-based taxonomic classifier in the spirit of Phymm.
    The main difference with Phymm is that we build IMMs for higher-level clades by pulling sequence
    data for subnodes.
    This class can be mostly viewed as imposing a TaxaTree structure onto ImmApp."""

    batchDepModes = ("predict","score","train","make-ref-seqdb",
            "fin-ref-seqdb","make-bench","bench-one-frag-len",
            "bench")

    ## Special taxonomy ID value to mark rejected samples 
    rejTaxid = rejTaxid #from MGT.TaxaConst
    
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
            --inp-train-seq 'test_data/seqdb-fasta/*.fasta.gz' \\
            --inp-train-model-descr 'test_data/seqdb-fasta/models.json' \\
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
                    "mode is used, the default will be the central database "+\
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
            
            optParseMakeOption_Path(None, "--imm-id-to-seq",
            dest="immIdToSeqIds",
            help="File that maps IMM IDs to lists of seq IDs during IMM training"),
            
            optParseMakeOption_Path(None, "--imm-id-to-meta",
            dest="immIdToMeta",
            help="File that maps IMM IDs to model meta data during IMM training"),
            
            optParseMakeOption_Path(None, "--imm-ids",
            dest="immIds",
            help="File with list of IMM IDs to use in scoring and prediction. Default is all"+\
                    " IMMs from --imm-seq-ids"),
            
            make_option(None, "--score-taxids-exclude-trees",
            action="store", 
            type="string",
            dest="scoreTaxidsExcludeTrees",
            help="String with a comma-separated taxonomy IDs of sub-trees to "+\
                    "exclude during scoring. This is never applied to --db-imm-archive collections."),
            
            optParseMakeOption_Path(None, "--inp-seq",
            dest="inpSeq",
            help="File or SeqDbFasta store with input FASTA sequences for prediction"),
            
            optParseMakeOption_Path(None, "--inp-seq-attrib",
            dest="sampAttrib",
            help="Optional tab-delimited file with extra attributes for each input sequence. "+\
                    "Currently must have two columns (no header row): sample id and weight. "+\
                    "Weight can be read count in a contig, and will be used when calculating "+\
                    "summary abundance tables."),

            optParseMakeOption_Path(None, "--inp-train-seq",
            dest="inpTrainSeq",
            help="File or shell glob with input FASTA sequences for training models"),
            
            optParseMakeOption_Path(None, "--inp-train-seq-list",
            dest="inpTrainSeqList",
            help="File that contains a list of file names (one per line) of the "+\
                    "input FASTA sequences for training models. It can also contain "+\
                    "a second tab-separated field with a shell glob pattern. "+\
                    "In that case, the first field can be either PathHasher directory "+\
                    "(see MGT/Util.py) or part of a file name preceding the glob."),
            
            make_option(None, "--inp-train-seq-format",
            action="store",
            type="choice",
            choices=("generic","ncbi"),
            default="generic",
            dest="inpTrainSeqFormat",
            help="Format of input training sequences: {ncbi,generic} [%default]. " +\
                    "'generic' requires --inp-train-model-descr option defined"),
            
            optParseMakeOption_Path(None, "--inp-train-model-descr",
            dest="inpTrainModelDescr",
            help="File in JSON format that maps models to training sequences"),
            
            make_option(None, "--inp-ncbi-seq-sel-policy",
            action="store",
            type="choice",
            choices=("drop-plasmids","extra-chrom-only","drop-extra-chrom"),
            default="drop-extra-chrom",
            dest="inpNcbiSeqSelPolicy",
            help="Policy preset for filtering NCBI training sequence based on the type of "+\
                    "genomic element [%default]"),
            
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
        
            optParseMakeOption_Path(None, "--out-id-score-meta",
            dest="outIdScoreMeta",
            help="Output file for combined meta-data map of score IDs "+\
            "[--out-score-comb+id_score.meta.pkl]"),
            
            optParseMakeOption_Path(None, "--taxa-tree-ncbi-dir",
            dest="taxaTreeNcbiDir",
            help="Directory where standard NCBI taxonomy tree is saved as a set of NCBI 'dump' files. "+\
                    "Mutually exclusive with --taxa-tree-pkl."),
            
            optParseMakeOption_Path(None, "--taxa-tree-pkl",
            dest="taxaTreePkl",
            help="Custom taxonomy tree saved in JSON format. If not set, standard NCBI tree is used, "+\
                    "except when training custom IMMs when this has to be always defined."),
            
            make_option(None, "--pred-min-len-samp",
            action="store", 
            type="int",
            dest="predMinLenSamp",
            help="Min length of samples to consider for prediction. 300 is default "+\
                    "for bacterial classification; 5000 is default for predicting "+\
                    "hosts for viral contigs."),
            
            optParseMakeOption_Path(None, "--pred-score-rescale-model",
            dest="predScoreRescaleModel",
            help="JSON file that describes the re-scaling model for scores "+\
                    "[data from global options]. This will be used to translate "+\
                    "raw prediction scores into normalized confidence scores."),
            
            make_option(None, "--train-composite-models",
            action="store", 
            type="int",
            default=0,
            dest="trainCompositeModels",
            help="If non-zero, train models for the inner nodes of the taxonomic "+\
                    "tree by concatenating balanced-picked leaf sequences [%default]."),
            
            make_option(None, "--train-min-len-samp",
            action="store", 
            type="int",
            dest="trainMinLenSamp",
            help="Min length of leaf node sequence to consider the node for training models: "+\
                    "default is 2000 if training custom models and 2000 if training "+\
                    "reference models. For training custom models, this means that sequences"+\
                    "shorter that this will be ignored"),
            
            make_option(None, "--train-min-len-samp-model",
            action="store", 
            type="int",
            dest="trainMinLenSampModel",
            help="Min length of sequence to use when training one model, otherwise the model "+\
                    "is skipped [10000]"),
            
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
            
            optParseMakeOption_Path(None, "--pred-out-stats-html",
            dest="predOutStatsHtml",
            help="Output HTML (Krona) file with statistics on predicted taxa [--pred-out-taxa/stats/stats.html]"),
            
            make_option(None, "--pred-out-stats-krona-url",
            action="store", 
            type="string",
            dest="predOutStatsKronaUrl",
            help="URL to Krona sources as seen by generated HTML output. See also --pred-out-stats-krona-embed"),
            
            make_option(None, "--pred-out-stats-krona-embed",
            action="store", 
            type="string",
            dest="predOutStatsKronaEmbed",
            help="Embed Krona sources at a given location. A relative location is interpreted " + \
                    "as being relative to --pred-out-stats-html. If --pred-out-stats-krona-url is not given, " + \
                    "its value is also computed from this option."),
            
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
            choices=("host","taxa","taxa-vir"),
            default="taxa",
            dest="predMode",
            help="Set the prediction mode [%default]: 'host' will work in a mode that assigns "+\
                    "host taxonomy to (presumed) bacteriophage "+\
                    "sequence or viral taxonomy; 'taxa' will assume a mixed sample "+\
                    "and assign taxonomy using all models; 'taxa-vir' will assume a "+\
                    "viral sample and use only viral models to make predictions."),
            
            make_option(None, "--skip-pred-out-taxa-csv",
            action="store",
            type="int",
            default=0,
            dest="skipPredOutTaxaCsv",
            help="If set to non-zero, --pred-out-taxa-csv file will not be produced - "+\
                    "useful when you only want the summary stats"),
            
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
            default="100,400,1000,5000,10000",
            dest="benchFragLenList",
            help="List of fragment lengths to generate and score in benchmark"),
            
            optParseMakeOption_Path(None, "--db-bench",
            dest="dbBench",
            help="BenchDb path, used either to create a new BenchDb from SeqDb "+\
                "or to use an existing BenchDb for benchmarking the models [-cwd/dbBench]"),
            
            optParseMakeOption_Path(None, "--db-bench-taxids-include",
            dest="dbBenchTaxidsInclude",
            help="File with taxonomy IDs, one per line, to "+\
                    "exclusively include during benchmark testing set construction."),
            
            make_option(None, "--db-bench-taxids-exclude-trees",
            action="store", 
            type="string",
            dest="dbBenchTaxidsExcludeTrees",
            help="String with a comma-separated taxonomy IDs of sub-trees to "+\
                    "exclude during benchmark testing set construction."),
            
            make_option(None, "--db-bench-filter-by-models",
            action="store",
            type="int",
            default=1,
            dest="dbBenchFilterByModels",
            help="If 1 [default], benchmark will only include samples "+\
                    "for taxa that have models in any of --db-imm. "+\
                    "Typically, this has to be set to 1 if you "+\
                    "are not using --bench-id-seq-db-to-id-score-remapping "+\
                    "option during processing of scoring results."),
            
            make_option(None, "--db-bench-frag-len",
            action="store",
            type="int",
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
            default=1,
            dest="benchNLevTestModelsMin",
            help="Minimum number of models below a node at a given taxonomic level required to evaluate "+\
                    "a benchmark performance for this node after all models related to the testing "+\
                    "sample were excluded at a given taxonomic exclusion level. The lowest possible "+\
                    "number is one, because otherwise we have no chance to produce the true prediction."),
            
            optParseMakeOption_Path(None, "--bench-taxids-map",
            dest="benchTaxidsMap",
            help="Depricated. Tab-delimited file with pairs input taxid, output taxid "+\
                    "that will be used to re-map taxonomy of the testing samples "+\
                    "before computing performance metrics. Used for phage-host benchmarking."),
            
            optParseMakeOption_Path(None, "--bench-id-seq-db-to-id-score-remapping",
            dest="benchIdSeqDbToIdScoreRemapping",
            help="JSON file that maps benchmark SeqDbIds to the known truth score IDs "+\
            "of the testing samples before computing performance metrics. "+\
            "Used for phage-host benchmarking."),
            
            optParseMakeOption_Path(None, "--bench-out-dir",
            dest="benchOutDir",
            help="Output directory for benchmarking results [-cwd/benchResults]"),
            
            optParseMakeOption_Path(None, "--bench-out-db-sqlite",
            dest="benchOutDbSqlite",
            help="Output SQLite database file with benchmarking results [--bench-out-dir/bench.sqlite]"),
            
            make_option(None, "--bench-proc-subtrees",
            action="append", 
            type="string",
            dest="benchProcSubtrees",
            help="Path of a file that contains taxonomy IDs, one per line. "+\
                "If set, benchmark performance metrics will be computed for "+\
                "subtrees of these IDs (inclusive) only."),
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
            if opt.mode in ("predict","make-bench","bench","proc-bench-scores"):
                opt.immDb = globOpt.icm.icmDbs
        if isinstance(opt.benchFragLenList,str):
            opt.benchFragLenList = [ int(x) for x in opt.benchFragLenList.split(",") ]
        #This options is accepted from the Web as a command line string, parsing should be done securily:
        if isinstance(opt.scoreTaxidsExcludeTrees,str):
            opt.scoreTaxidsExcludeTrees = [ int(x) for x in opt.scoreTaxidsExcludeTrees.split(",") ]
            #process special case for Web interface
            if opt.scoreTaxidsExcludeTrees == [0]:
                opt.scoreTaxidsExcludeTrees = None
        optPathMultiOptToAbs(opt,"immDb")           
        optPathMultiOptToAbs(opt,"immDbArchive")           
        opt.setIfUndef("immIdToSeqIds",pjoin(opt.cwd,"imm-id-to-seq"))
        opt.setIfUndef("immIdToMeta",pjoin(opt.cwd,"imm-id-to-meta"))
        opt.setIfUndef("outDir",pjoin(opt.cwd,"scores"))
        opt.setIfUndef("predOutDir",pjoin(opt.cwd,"results"))
        opt.setIfUndef("seqDb",pjoin(opt.cwd,"seqDb"))
        opt.setIfUndef("dbBench",pjoin(opt.cwd,"dbBench"))
        opt.setIfUndef("immIds",opt.immIdToSeqIds)
        opt.setIfUndef("outScoreComb",pjoin(opt.outDir,"combined"+ImmApp.scoreSfx))
        opt.setIfUndef("outIdScoreMeta",opt.outScoreComb+".id_score.meta.pkl")
        opt.setIfUndef("predOutTaxa",pjoin(opt.predOutDir,"pred-taxa"))
        opt.setIfUndef("predOutTaxaCsv",opt.predOutTaxa+".csv")
        opt.setIfUndef("predOutTaxaSummCsv",opt.predOutTaxa+".summ.csv")
        opt.setIfUndef("predOutDbSqlite",opt.predOutTaxa+".sqlite")
        opt.setIfUndef("writeEachSampleSqlite",False)
        opt.setIfUndef("predOutStatsDir", pjoin(opt.predOutDir,"stats"))
        opt.setIfUndef("predOutStatsCsv", pjoin(opt.predOutStatsDir,"stats.csv"))
        opt.setIfUndef("predOutStatsPdf", pjoin(opt.predOutStatsDir,"stats.pdf"))
        opt.setIfUndef("predOutStatsHtml", pjoin(opt.predOutStatsDir,"stats.html"))
        opt.setIfUndef("newTaxidTop",mgtTaxidFirst)
        opt.setIfUndef("immDbWorkDir",pjoin(opt.cwd,"immDbWorkDir"))
        opt.setIfUndef("scoreWorkDir",pjoin(opt.cwd,"scoreWorkDir"))
        opt.setIfUndef("benchWorkDir",pjoin(opt.cwd,"benchWorkDir"))
        opt.setIfUndef("benchOutDir",pjoin(opt.cwd,"benchResults"))
        opt.setIfUndef("benchOutDbSqlite",pjoin(opt.benchOutDir,"bench.sqlite"))
        optPathMultiOptToAbs(opt,"benchProcSubtrees")

        if not opt.isUndef("predOutStatsKronaEmbed"):
            opt.setIfUndef("predOutStatsKronaUrl",opt.predOutStatsKronaEmbed)


        if isinstance(opt.dbBenchTaxidsExcludeTrees,str):
            opt.dbBenchTaxidsExcludeTrees = [ int(x) for x in opt.dbBenchTaxidsExcludeTrees.split(",") ]
        
        opt.setIfUndef("rejectRanksHigher","superkingdom")
        
        if opt.predMode == "host": 
            opt.setIfUndef("predMinLenSamp",5000)
        elif opt.predMode == "taxa":
            opt.setIfUndef("predMinLenSamp",300)
        elif opt.predMode == "taxa-vir":
            opt.setIfUndef("predMinLenSamp",300)
            opt.setIfUndef("scoreTaxidsExcludeTrees",[cellTaxid])
        else:
            raise ValueError("Unknown --pred-mode value: %s" % (opt.predMode,))
        
        if opt.mode == "train":
            opt.setIfUndef("trainMinLenSamp",2000)
            opt.setIfUndef("trainMinLenSampModel",10000)

        if opt.mode == "make-ref-seqdb":
            globOpt = globals()["options"]
            if opt.isUndef("inpTrainSeq") and opt.isUndef("inpTrainSeqList"):
                opt.inpTrainSeqFormat = "ncbi"
                ph = PathHasher(globOpt.refSeqDataDir,mode="r")
                opt.inpTrainSeqList = pjoin(opt.cwd,"inp_train_seq.csv")
                with open(opt.inpTrainSeqList,"w") as out:
                    for f in ph.glob("*.fna.gz"):
                        out.write("{}\n".format(f))
            if opt.inpTrainSeqFormat == "generic":
                if opt.isUndef("inpTrainModelDescr"):
                    raise ValueError("--inp-train-model-descr must be defined when --inp-train-seq-format is 'generic'")
        if opt.benchNLevTestModelsMin < 1:
            parser.error("--bench-n-lev-test-model-min must be at least 1")

    
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
        
        #these levels will not be reported in the prediction output -
        #"root" is trivial and "kingdom" - for consistency between 
        #euks and proks
        self.skipLevNamesInPredOut = ("kingdom","root")
    
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
            return self.setupTrainingSimple(**kw)
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
            self.taxaTree = loadTaxaTree(ncbiTaxaDumpDir=self.opt.taxaTreeNcbiDir,
                    jsonFile=self.opt.taxaTreePkl)
        return self.taxaTree

    def getTaxaLevels(self):
        if self.taxaLevels is None:
            #that assigns "level" and "idlevel" attributes to TaxaTree nodes,
            #when taxaTree is already loaded. Otherwise, you can use
            #taxaLevels.setTaxaTree() later.
            self.taxaLevels = TaxaLevels(self.taxaTree)
        return self.taxaLevels
    
    def getSeqDb(self,mode="r"):
        opt = self.opt
        return SeqDbFasta.open(opt.seqDb,mode=mode)

    def makeRefSeqDb(self,**kw):
        """Create reference SeqDb (the "main" SeqDb)"""
        opt = self.opt
        optI = copy(opt)
        optI.mode = "extract-ref-seqdb"
        #We need to finish the extraction now, because
        #the fin-ref-seqdb stage needs to know all IDs 
        #collected. Extraction is not parralelized anyway.
        optI.runMode = "inproc"
        app = self.factory(opt=optI)
        jobs = app.run(**kw)
        optI = copy(opt)
        optI.mode = "fin-ref-seqdb"
        app = self.factory(opt=optI)
        kwI = kw.copy()
        if jobs:
            kwI["depend"] = jobs
        return app.run(**kwI)
    
    def extractRefSeqDb(self,**kw):
        """Extract reference SeqDb from input database sequences"""
        opt = self.opt
        self.seqDb = None
        seqDb = SeqDbFasta.open(path=opt.seqDb,mode="c")
        taxaTree = self.getTaxaTree()
        seqDb.setTaxaTree(taxaTree)
        inpTrainSeqFiles = []
        if not opt.isUndef("inpTrainSeq"):
            #these can be only regular files or patterns by not PathHasher
            #because we would know the glob (e.g. extension) inside PathHasher
            inpTrainSeqFiles += list(glob.glob(os.path.abspath(opt.inpTrainSeq)))
        if not opt.isUndef("inpTrainSeqList"):
            with closing(openCompressed(opt.inpTrainSeqList,"r")) as inp:
                for line in inp:
                    line = line.strip()
                    if line:
                        words = line.split("\t")
                        rootPath = os.path.abspath(words[0])
                        if len(words) == 1:
                            #just a file name
                            inpTrainSeqFiles += rootPath
                        else:
                            #root path\tpattern
                            patt = words[1]
                            if PathHasher.is_instance_dir(rootPath):
                                #root path is PathHasher
                                ph = PathHasher(rootPath,mode="r")
                                inpTrainSeqFiles += [ f for f in ph.glob(patt) ]
                            else:
                                #root path with pattern is just a pattern
                                fullPatt = pjoin(rootPath,patt)
                                inpTrainSeqFiles += list(glob.glob(fullPatt))
       
        if opt.inpTrainSeqFormat == "ncbi":
            policyFilter=FastaTrainingSeqFilter(policy=opt.inpNcbiSeqSelPolicy)
            filt = functools.partial(fastaReaderFilterNucDegen,
                    extraFilter=policyFilter,
                    minNonDegenRatio=0.95)
            splitFastaFilesByTaxa(inSeqs=inpTrainSeqFiles,
                    outStore=seqDb,
                    taxaTree=taxaTree,
                    filt=filt)

        elif opt.inpTrainSeqFormat == "generic":
            filt = functools.partial(fastaReaderFilterNucDegen,
                    minNonDegenRatio=0.90)
            splitFastaFilesByModel(inSeqs=inpTrainSeqFiles,
                    modelsMeta=loadTrainModelsDescr(opt.inpTrainModelDescr),
                    outStore=seqDb,
                    taxaTree=taxaTree,
                    checkTaxa=True,
                    filt=filt)


        # For now, we filter by length when we build models.
        #taxids = seqDb.getIdList()
        #for taxid in taxids:
        #    taxidSeqLen = seqDb.seqLengths(taxid)["len"].sum()
        #    if taxidSeqLen < opt.trainMinLenSamp:
        #        seqDb.delById(taxid)

    def finRefSeqDb(self,**kw):
        """Finalize creation of SeqDb for training ICM models"""
        opt = self.opt
        seqDb = self.getSeqDb(mode="r+")
        ids = seqDb.getIdList(objSfx=seqDb.objUncomprSfx)
        ids = n.asarray(ids,dtype="O")
        nrnd.shuffle(ids)
        jobs = []
        for idsBatch in n.array_split(ids,min(1000,len(ids))):
            optI = copy(opt)
            optI.mode = "fin-ref-seqdb-batch"
            optI.seqDbIds = idsBatch
            app = self.factory(opt=optI)
            jobs += app.run(**kw)
        return jobs
    
    def finRefSeqDbBatch(self,**kw):
        """Sub-task of finalizing creation of SeqDb for training ICM models"""
        opt = self.opt
        seqDb = self.getSeqDb(mode="w")
        ids = opt.seqDbIds
        for id in ids:
            seqDb.finById(id)
            #print "DEBUG: finByid(%s) done" % (taxid,)

    ## Methods that generate benchmark dataset

    def makeBenchFromSeqDb(self,**kw):
        """Generate a benchmark dataset from SeqDb that was used to train ICM models"""
        opt = self.opt
        seqDb = self.getSeqDb()
        bench = ImmClassifierBenchmark.open(path=opt.dbBench,mode="c")
        bench.seqDb = seqDb
        if opt.dbBenchFilterByModels:
            immDbs = [ ImmStoreWithTaxids(immDb) for immDb in opt.immDb ]
        else:
            immDbs = None
        ids = bench.selectIdsDb(immDbs=immDbs,
                dbBenchTaxidsExcludeTrees=opt.dbBenchTaxidsExcludeTrees,
                dbBenchTaxidsInclude=opt.dbBenchTaxidsInclude)
        ids = n.asarray(list(ids),dtype="O")
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
        bench = ImmClassifierBenchmark.open(path=opt.dbBench,mode="r+")
        bench.seqDb = seqDb
        for id in ids:
            bench.makeSample(idDb=id,
                    fragLen=opt.dbBenchFragLen,
                    fragCountMax=opt.dbBenchFragCountMax)
        
    def makeBenchFromSeqDbFin(self,**kw):
        """Final task of benchmark generation"""
        return kw.get("depend",None)

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
            #see comment below for why we creating an empty opt object
            optI = Struct()
            optI.runMode = opt.runMode
            optI.lrmUserOptions = opt.lrmUserOptions
            optI.mode = "bench-one-frag-len"
            optI.immDb = opt.immDb
            optI.seqDb = opt.seqDb
            optI.dbBenchFragLen = fragLen
            optI.dbBenchFragCountMax = opt.dbBenchFragCountMax
            optI.dbBenchFilterByModels = opt.dbBenchFilterByModels
            optI.dbBenchTaxidsExcludeTrees = opt.dbBenchTaxidsExcludeTrees
            optI.dbBenchTaxidsInclude = opt.dbBenchTaxidsInclude
            optI.scoreTaxidsExcludeTrees = opt.scoreTaxidsExcludeTrees
            optI.benchNLevTestModelsMin = opt.benchNLevTestModelsMin
            optI.benchTaxidsMap = opt.benchTaxidsMap
            optI.benchIdSeqDbToIdScoreRemapping = opt.benchIdSeqDbToIdScoreRemapping
            optI.benchProcSubtrees = opt.benchProcSubtrees
            optI.cwd = pjoin(opt.benchWorkDir,str(fragLen))
            optI.benchOutDir = pjoin(opt.benchOutDir,str(fragLen))
            optI.dbBench = pjoin(opt.dbBench,str(fragLen))
            optI.nImmBatches = opt.nImmBatches
            #this sets default and derived options
            #@todo we need to create two kinds of properties in Struct():
            #fixed values and lazy-evaluated values (similar to variables
            #in "make" assigned with "=" and ":="). This will allow changing
            #just one option in a copy of the parent option object (e.g.
            #the working directory), with dependent options (e.g. various
            #filenames) dynamically recomputed.
            #This is currently a reason that we cannot just copy the opt above - 
            #we would get all absolutized output paths set to identical values
            #for each of the frag-len jobs.
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
        optI.inpSeq = opt.dbBench
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

    def _makeTaxaSampLengths(self,seqDb,trainMinLenSamp):
        taxaSampLengths = {} 
        for seqid in seqDb.getIdList():
             taxaSampLen = seqDb.seqLengthSum(seqid) 
             if taxaSampLen >= trainMinLenSamp:
                 taxaSampLengths[seqid] = taxaSampLen
        return taxaSampLengths

    def mapSeqToTree(self):
        """Assign list of SeqDB IDs to corresponding taxonomy tree nodes.
        In the current SeqDB version, each ID is a unique taxonomy id, so
        a one-element list will be assigned.
        @post Attribute 'leafSeqDbIds' is assigned to EVERY node and contains a list of IDs, possibly empty.
        The empty list will be a reference to a single shared object, so it should be treated as immutable"""
        opt = self.opt
        taxaTree = self.getTaxaTree()
        seqDb = self.getSeqDb()
        taxaSampLengths = self._makeTaxaSampLengths(seqDb=seqDb,trainMinLenSamp=opt.trainMinLenSamp)
        self.taxaSampLengths = taxaSampLengths
        emptyList = []
        taxaTree.setAttribute("leafSeqDbIds",emptyList,doCopy=False)
        for taxid in taxaSampLengths:
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
                    #at the bottom of the lineage. This implicitly selects the
                    #most specific assignment for such lineages, unless we change
                    #it at the prediction stage by looking at the location of models
                    #on the tree.
                    if len(node.leafSeqDbIds)>0:
                        node.trainSelStatus = self.TRAIN_SEL_STATUS_DIRECT
                    elif len(chSeqIds) >= 2:
                        node.trainSelStatus = self.TRAIN_SEL_STATUS_INDIRECT


        taxaTree.visitDepthBottom(actor)

    def defineImms(self):
        opt = self.opt
        taxaTree = self.getTaxaTree()
        micNodes = [ taxaTree.getNode(taxid) for taxid in micTaxids ]
        cellNode = taxaTree.getNode(cellTaxid)
        virNode = taxaTree.getNode(virTaxid)
        micVirNodes = micNodes+[virNode]
        #DEBUG:
        #cntLeafSeq = sum([ len(node.leafSeqDbIds)>0 for node in taxaTree.iterDepthTop() if node.isUnder(cellNode) ])
        #cntDirect = sum([ hasattr(node,"trainSelStatus") and node.trainSelStatus == self.TRAIN_SEL_STATUS_DIRECT for node in taxaTree.iterDepthTop() if node.isUnder(cellNode) ])
        immIdToSeqIds = {}
        taxaSampLengths = self.taxaSampLengths
        ##@todo Make it controlled by a node selection algebra passed by the user as json expression
        for node in taxaTree.iterDepthTop():
            if hasattr(node,"trainSelStatus") and node.trainSelStatus != self.TRAIN_SEL_STATUS_IGNORE:
                doPick = False
                if not opt.trainCompositeModels:
                    if node.trainSelStatus == self.TRAIN_SEL_STATUS_DIRECT:
                        doPick = True
                else:
                    if node.isUnderAny(micVirNodes):
                        doPick = True
                    elif node.isUnder(cellNode):
                        if node.trainSelStatus == self.TRAIN_SEL_STATUS_DIRECT:
                            doPick = True
                if doPick:
                    nodeSampLength = sum(taxaSampLengths[taxid] for taxid in node.pickedSeqDbIds)
                    if nodeSampLength >= opt.trainMinLenSampModel:
                        #DEBUG:
                        #print "Training for: ", node.lineageStr()
                        immIdToSeqIds[node.id] = node.pickedSeqDbIds
        dumpObj(immIdToSeqIds,opt.immIdToSeqIds)
        immIdToMeta = dict( 
                ( (taxid,
                    dict(taxid = taxid,
                        is_leaf = len(seqIds) == 1,
                        name = taxaTree.getNode(taxid).name
                        )
                    ) for (taxid,seqIds) in immIdToSeqIds.items()) )
                #TODO: no longer create new nodes; use this function only for DB sequences; or, rather, assign idModel[] to nodes and use those rather than node itself
        dumpObj(immIdToMeta,opt.immIdToMeta)

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
        Each scaffold is treated as a separate taxonomic unit under the super-node of environmental bacterial sequences.
        Alternatively, if the FASTA defline starts with TAXIDXXX_, where XXX 
        is an integer number, the number is interpreted as a taxid of the
        parent node. In that case the parent taxid must be already present in
        the currently loaded tree.
        This method recreates the currently loaded tree."""
        opt = self.opt
        assert opt.taxaTreePkl, "File name for a custom taxonomic tree must be provided"
        self.seqDb = None
        seqDb = SeqDbFasta.open(path=opt.seqDb,mode="c")
        compr = SymbolRunsCompressor(sym="Nn",minLen=1)
        nonDegenSymb = "ATCGatcg"
        newTaxaTop = TaxaNode(id=opt.newTaxidTop,name=opt.newTaxNameTop,
                rank=unclassRank,divid=dividEnv,names=list())
        newTaxaTopUsed = False
        nextNewTaxid = newTaxaTop.id + 1
        fastaReader = FastaReader(opt.inpTrainSeq)
        nNodesOut = 0
        #load pristine NCBI tree because we will save it at the end
        #after adding new nodes
        self.taxaTree = None
        taxaTree = loadTaxaTree(ncbiTaxaDumpDir=opt.taxaTreeNcbiDir)
        self.taxaTree = taxaTree
        #####
        for rec in fastaReader.records():
            hdr = rec.header()
            seqid = rec.getId() # scf768709870 or TAXIDXXXX_YYYYY
            seq = compr(rec.sequence())
            if not checkSaneAlphaHist(seq,nonDegenSymb,minNonDegenRatio=0.90):
                print "WARNING: ratio of degenerate symbols is too high, "+\
                        "skipping the reference scaffold id %s" % (seqid,)
            if len(seq) >= opt.trainMinLenSamp:
                taxaNode = TaxaNode(id=nextNewTaxid,name=self._customTrainSeqIdToTaxaName(seqid),rank=unclassRank,divid=dividEnv,names=list())
                if seqid.lower().startswith("taxid"):
                    parTaxid = int(seqid[len("taxid"):].split("_",1)[0])
                    try:
                        parent = taxaTree.getNode(parTaxid)
                    except:
                        print "Parent taxid {0} not found".format(seqid)
                        parent = taxaTree.getRootNode()
                else:
                    parent = newTaxaTop
                    newTaxaTopUsed = True
                taxaNode.setParent(parent)
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
        taxaSampLengths = self._makeTaxaSampLengths(seqDb=seqDb,trainMinLenSamp=opt.trainMinLenSamp)
        self.taxaSampLengths = taxaSampLengths
        if newTaxaTopUsed:
            envBacTop = taxaTree.getNode(envBacTaxid)
            newTaxaTop.setParent(envBacTop)
        taxaTree.rebuild()
        taxaTreeStore = NodeStorageJson(opt.taxaTreePkl)
        taxaTreeStore.save(taxaTree)
    
    def defineImmsSimple(self,**kw):
        opt = self.opt
        taxaTree = self.getTaxaTree()
        seqDb = self.getSeqDb()
        immIdToSeqIds = dict()
        immIdToMeta = dict()
        for idSeqDb in seqDb.getIdList():
            sampLen = seqDb.seqLengthSum(idSeqDb)
            if sampLen >= min(opt.trainMinLenSamp,opt.trainMinLenSampModel):
                seqMeta = seqDb.loadMetaDataById(idSeqDb)
                modMeta = dict(
                    taxid = seqMeta["taxid"],
                    name = seqMeta["name"],
                    is_leaf = 1)
                idMod = idSeqDb
                immIdToMeta[idMod] = modMeta
                immIdToSeqIds[idMod] = [idSeqDb]
        dumpObj(immIdToSeqIds,opt.immIdToSeqIds)
        dumpObj(immIdToMeta,opt.immIdToMeta)
    
    def setupTrainingSimple(self,**kw):
        return self.defineImmsSimple(**kw)

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
            optI.path = dbArch[0]
            optI.archive = dbArch[1]
            if opt.web:
                app = ArchiveGpgApp(opt=optI)
            else:
                optI.stripComponents = 1
                app = ArchiveApp(opt=optI)
            kw = kw.copy()
            kw["depend"] = jobs
            jobs = app.run(**kw)
        
        return jobs

    def score(self,**kw):
        """Score with all IMMs.
        Parameters are taken from self.opt
        @param inpSeq Name of the input multi-FASTA file or SeqDbFasta store to score
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        @param outIdScoreMeta name for output file with combined meta
        data for each idScore (idScore equals idImm in this class)
        """
        opt = self.opt

        optI = copy(opt)
        optI.mode = "split"
        optI.chunkSize = 100*1024**2
        scoreSeqRoot = pjoin(opt.cwd,"inp-seq-scoring.db")
        optI.outSeq = scoreSeqRoot
        optI.outSeqDbIds = scoreSeqRoot+".ids"
        optI.assertOutSeq = 1
        app = SeqDbFastaApp(opt=optI)
        kwI = kw.copy()
        jobsD = app.run(**kwI)

        opt.inpSeq = optI.outSeq
        opt.inpSeqDbIds = optI.outSeqDbIds

        immDb = [ (d,None) for d in opt.immDb ]
        immDbArch = self._archiveNamesToDirNames(opt.immDbArchive,opt.immDbWorkDir,"imm")
        immDb += immDbArch
        outSubScores = []
        idScoreMeta = dict()
        jobs = []
        for immD in immDb:
            jobsI = copy(jobsD)
            if immD[1] is not None:
                optI = copy(opt)
                optI.mode = "extract"
                #we need to get metadata, so need to expand here and now
                #@todo maybe extract only metadata dir, and then submit
                #full extraction to batch
                optI.runMode = "inproc"
                optI.path = immD[0]
                optI.archive = immD[1]
                if opt.web:
                    optI.tarArgs = "--strip-components 1"
                    app = ArchiveGpgApp(opt=optI)
                else:
                    app = ArchiveApp(opt=optI)
                kwI = kw.copy()
                #jobsI = app.run(**kwI)
                #immIds = ImmStoreWithTaxids(immD[0]).listImmIdsWithTaxids(
                #        iterPaths=(item.name for item in app.iterMembers()))
                app.run(**kwI)
                immIds = ImmStoreWithTaxids(immD[0]).dictMetaData()
                assert len(immIds) > 0,"No IMMs found in IMM DB - probably training did not run yet"
            else:
                immIds = ImmStoreWithTaxids(immD[0]).dictMetaData()
                assert len(immIds) > 0,"No IMMs found in IMM DB - probably training did not run yet"
                #we only aply taxid filters to "reference" immDbs, under 
                #an assumption that archived immDbs are supplied by the user
                #and contain only what needs to be used
                if opt.scoreTaxidsExcludeTrees is not None:
                    taxaTree = self.getTaxaTree()
                    try:
                        subTreesExcl = [ taxaTree.getNode(taxid) for taxid \
                                in set(opt.scoreTaxidsExcludeTrees) ]
                    except KeyError,msg:
                        raise ValueError,"Unknown taxid: %s" % (msg,)
                    immIds = dict(( (id,meta) for (id,meta) in immIds.items() \
                            if not taxaTree.getNode(meta["taxid"]).isUnderAny(subTreesExcl) )) 
                    if len(immIds) <= 0:
                        print "Warning: No IMMs left in IMM DB after applying exclusion filters"

            
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
            #define idScore == idImm:
            dumpObj([(immId,immId) for immId in immIds],optI.immIds)
            optI.outScoreComb = pjoin(optI.outDir,"combined"+ImmApp.scoreSfx)
            outSubScores.append(optI.outScoreComb)
            #define idScore == idImm:
            idScoreMeta.update(immIds)
            app = ImmApp(opt=optI)
            jobsI = app.run(**kwI)
            jobs += jobsI
        
        makeFilePath(opt.outIdScoreMeta)
        dumpObj(idScoreMeta,opt.outIdScoreMeta)

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
        outScoreCombWork = makeWorkFile(opt.outScoreComb)
        immScores = openImmScores(opt,fileName=outScoreCombWork,mode="w")
        immScores.catScores(fileNames=opt.outSubScores)
        immScores.close()
        os.rename(outScoreCombWork,opt.outScoreComb)

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

    def _roundUpPredictions(self,predTaxid,rootNodes,rank):
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
        rankPos = taxaLevels.getLevelPos(rank)
        for (i,taxid) in enumerate(predTaxid):
            if not taxid == self.rejTaxid:
                node = taxaTree.getNode(taxid)
                if node.isSubnodeAny(rootNodes):
                    lin = taxaLevels.lineageFixedList(node,null=None,format="node",fill="up-down")
                    rankNode = lin[rankPos]
                    #assert rankNode is not None,("Logic error for predNode=%s" % (predNode,))
                    #condition below should probably always be true...
                    if rankNode:
                        predTaxid[i] = rankNode.id

    def scoresToPredictions(self,immScores,idScoreMeta,**kw):
        """Generate taxonomic predictions from raw IMM scores.
        This handles both classification of bacterial sequences and host assignment 
        for viral sequences.
        @param immScores ImmScores object.
        @param idScoreMeta dict(idScore->metadata)
        
        Other parameters are taken from self.opt.
        @param rndScoreComb File with ImmScores object for random query sequences
        @param predOutTaxa Output file with predicted taxa
        """
        opt = self.opt
        sc = immScores.getData()
        makeFilePath(opt.predOutTaxa)
        pred = TaxaPred(opt.predOutTaxa,mode="w")
        pred.initFromSamp(idSamp=sc.idSamp,lenSamp=sc.lenSamp,idScoreDtype=sc.idScore.dtype)
        pr = pred.getData()
        _idScoreToTaxid = lambda id_scores: \
                [ idScoreMeta[id_score]["taxid"] for id_score in id_scores ]
        kind = immScores.getKind()
        if kind == "ImmScoresReduced":
            #implicit conversion into int
            array_chunked_copy(sc.idScore,pr.predTaxid,chunk=10**6,op=_idScoreToTaxid)
            array_chunked_copy(sc.idScore,pr.predIdScore,chunk=10**6)
            array_chunked_copy(sc.score,pr.predScore,chunk=10**6)
        elif kind == "ImmScoresDenseMatrix":
            #scRnd = loadObj(opt.rndScoreComb)
            #sc.score = self._normalizeScores(sc,scRnd)
            #normalize to Z-score along each row
            #sc.score = ((sc.score.T - sc.score.mean(1))/sc.score.std(1)).T
            #normalize to Z-score over entire matrix
            #sc.score = ((sc.score - sc.score.mean())/sc.score.std())
            #Until we optimize chunk layout of the HDF5 score matrix
            #(by using dataset copy method),
            #we have to load the entire matrix into RAM, otherwise
            #reading it by columns is extremely slow. Once this is
            #optimized, simply remove [:] after sc.score. Better
            #still, convert argmax to iteration to avoid creating
            #argmax array.
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
            #Implicit conversion to int
            predIdScore =  sc.idScore[argmaxSc]
            array_chunked_copy(predIdScore,
                    pr.predTaxid,chunk=10**6,op=_idScoreToTaxid)
            array_chunked_copy(predIdScore,
                    pr.predIdScore,chunk=10**6)
            array_chunked_copy(score[argmaxSc,n.arange(score.shape[1])],
                    pr.predScore,chunk=10**6)

        else:
            raise ValueError("Unknown kind of score object: %s" % (kind,))
        return pred

    def filterPredictions(self,immScores,predTaxid,**kw):
        """Generate taxonomic predictions from raw IMM scores.
        This handles both classification of bacterial sequences and host assignment 
        for viral sequences.
        @param immScores ImmScores object.
        @param predTaxid[in,out] array with predicted taxonomy per sample
        
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
        if opt.rejectRanksHigher and opt.rejectRanksHigher != "superkingdom":
            max_linn_levid = taxaLevels.getLevelId(opt.rejectRanksHigher)
            min_linn_levid = taxaLevels.getLinnLevelIdRange()[0]
            # Reject predictions to clades outside of certain clade level range
            for i in xrange(len(predTaxid)):
                if not predTaxid[i] == self.rejTaxid:
                    predNode = taxaTree.getNode(predTaxids[i])
                    # we need to protect leaf nodes because we place environmental scaffolds
                    # as no_rank under bacteria->environmental
                    if not (predNode.isLeaf() or taxaLevels.isNodeInLinnLevelRange(predNode,
                            min_linn_levid,max_linn_levid)):
                        predTaxid[i] = self.rejTaxid

        #Round-up euk predictions to the phylum level because of high sequence identities between
        #lower order clades
        self._roundUpPredictions(predTaxid,rootNodes=(taxaTree.getNode(eukTaxid),),rank="phylum")
    
    def processImmScores(self,**kw):
        """Process raw IMM scores to predict taxonomy.
        This handles both classification of bacterial sequences and host assignment 
        for viral sequences.
        Parameters are taken from self.opt.
        @param outScoreComb File with ImmScores object
        @param outIdScoreMeta File with meta data for each score ID
        @param rndScoreComb File with ImmScores object for random query sequences
        @param predOutTaxa Output file with predicted taxa
        """
        opt = self.opt
        immScores = openImmScores(opt,fileName=opt.outScoreComb,mode="r")
        idScoreMeta = loadObj(opt.outIdScoreMeta)
        pred = self.scoresToPredictions(immScores=immScores,
                idScoreMeta=idScoreMeta,
                **kw)
        self.filterPredictions(immScores=immScores,
                predTaxid=pred.getData().predTaxid,
                **kw)
        pred.close()
        self.exportPredictions(idScoreMeta=idScoreMeta,**kw)
        #todo: matrix of predhost idlevel's; batch imm.score() jobs; score hist;

    def exportPredictions(self,idScoreMeta,**kw):
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
        self._exportPredictions(idScoreMeta=idScoreMeta)
        #if you need all individual predictions as SQL too, set the option to True
        if opt.writeEachSampleSqlite:
            self._loadPredictionsIntoSql(csvFile=opt.predOutTaxaCsv,isSumm=False)
        self._loadPredictionsIntoSql(csvFile=opt.predOutTaxaSummCsv,isSumm=True)
        self.statsPred()

    def _buildCustomTaxidToRefTaxidMap(self,predOutTaxa):
        """Helper method for transitive classification.
        This loads a prediction file created in a separate
        run for assigning custom model sequences to a reference DB,
        and returns a dict that maps custom taxonomic id generated
        here to assigned reference DB taxid"""
        opt = self.opt
        pred = TaxaPred(predOutTaxa,mode="r")
        taxaPred = pred.getData()
        custToRef = {}
        for idS,taxidRef,lenS in it.izip(taxaPred.idSamp,taxaPred.predTaxid,taxaPred.lenSamp):
            taxidCust = self._customTrainSeqIdToTaxid(idS)
            if taxidCust:
                custToRef[taxidCust] = taxidRef
        pred.close()
        return custToRef

    def _exportPredictions(self,idScoreMeta):
        opt = self.opt    
        pred = TaxaPred(opt.predOutTaxa,mode="r")
        taxaPred = pred.getData()
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        levNames = taxaLevels.getLevelNames("ascend")
        skipLevNames = self.skipLevNamesInPredOut
        levIdRoot = taxaLevels.getLevelId("root")
        if opt.transPredOutTaxa:
            transTaxaMap = self._buildCustomTaxidToRefTaxidMap(opt.transPredOutTaxa)
        else:
            transTaxaMap = None
        if not opt.skipPredOutTaxaCsv:
            makeFilePath(opt.predOutTaxaCsv)
            out = openCompressed(opt.predOutTaxaCsv,"w")
        else:
            out = null_file()
        makeFilePath(opt.predOutStatsHtml)
        resc = ImmScoreRescaler(predScoreRescaleModel=opt.predScoreRescaleModel)
        krona = KronaWriter(taxaTree=taxaTree)
        krona.initCount()
        flds = ["id","len","taxid","name","rank","weight","idscore","namescore"]
        for levName in levNames:
            if levName not in skipLevNames:
                flds += ["taxid_"+levName,"name_"+levName,"score_"+levName]
        w = csv.DictWriter(out, fieldnames=flds, restval='null',dialect='excel-tab')
        w.writerow(dict([(fld,fld) for fld in flds]))
        def _row_iter_2d(arr,n_col):
            for row in arr:
                yield row[:n_col]
        def _row_iter_1d(arr):
            for el in arr:
                yield (el,)
        predPerSample = 1
        if predPerSample > 1:
            predTaxidRows = _row_iter_2d(taxaPred.topTaxid,predPerSample)
            predIdScoreRows = _row_iter_2d(taxaPred.topIdScore,predPerSample)
            predScoreRows = _row_iter_2d(taxaPred.topScores,predPerSample)
        else:
            predTaxidRows = _row_iter_1d(taxaPred.predTaxid)
            predIdScoreRows = _row_iter_1d(taxaPred.predIdScore)
            predScoreRows = _row_iter_1d(taxaPred.predScore)
        if opt.sampAttrib:
            sampAttribCsvFileInp = openCompressed(opt.sampAttrib,"r")
            sampAttribCsv = csv.reader(sampAttribCsvFileInp,dialect="excel-tab")
        else:
            sampAttribCsvFileInp = None
            sampAttribCsv = None
        idScoreSumm = {}
        predIter = it.izip(taxaPred.idSamp,predTaxidRows,
                                predScoreRows,taxaPred.lenSamp,
                                predIdScoreRows)
        if sampAttribCsv:
            recIter = join_sorted_records_right_check(iter1=predIter,
                    iter2=sampAttribCsv,
                    key1=0,
                    key2=0)
        else:
            recIter = it.izip(predIter,it.repeat(None))
        for ((idS,taxids,scores,lenS,idscores),recAttr) in recIter:
            if lenS < opt.predMinLenSamp:
                continue
            if recAttr:
                weight = float(recAttr[1])
            else:
                #if weight attribute is not provided, we use length of each sample
                #as weight
                weight = float(lenS)
            for (taxid,score,idscore) in it.izip(taxids,scores,idscores):
                row = dict(id=idS,len=lenS)
                if transTaxaMap:
                    if taxid != self.rejTaxid:
                        taxid = transTaxaMap[taxid]
                scoresResc = {}
                if taxid != self.rejTaxid:
                    node = taxaTree.getNode(taxid)
                    row["name"] = node.name
                    row["rank"] = node.linn_level
                    lin = taxaLevels.lineageFixedList(node,null=None,format="node",fill="up-down")
                    for (levName,levNode) in it.izip(levNames,lin):
                        if levNode:
                            if levName not in skipLevNames:
                                levIdReal = taxaLevels.getLevelId(levNode.linn_level)
                                assert levIdReal >= taxaLevels.minLinnId,\
                                        "Found norank node in supposedly Linnean lineage: %s" % (str(levNode),)
                                row["name_"+levName] = levNode.name
                                row["taxid_"+levName] = levNode.id
                                if levIdReal == levIdRoot:
                                    #this should not happen but just in case
                                    scoreResc = 1.
                                else:
                                    scoreResc = round(resc.rescale(score=score,
                                        len_samp=lenS,
                                        i_lev_pre=levIdReal),2)
                                scoresResc[levNode.id] = scoreResc
                                row["score_"+levName] = scoreResc
                row["taxid"] = taxid
                row["idscore"] = idscore
                row["namescore"] = idScoreMeta[idscore]["name"]
                row["weight"] = weight
                w.writerow(row)
                krona.addSample( (idS,taxid,weight,scoresResc) )
                #aggregate into taxa summary table
                if not idscore in idScoreSumm:
                    rec = dict(row)
                    rec["cnt"] = 1
                    idScoreSumm[idscore] = rec
                else:
                    rec = idScoreSumm[idscore]
                    rec["cnt"] += 1
                    rec["weight"] += row["weight"]
                    rec["len"] += row["len"]
                    for fld_name in rec.keys():
                        if fld_name.startswith("score_"):
                            rec[fld_name] += row[fld_name]
        out.close()
        pred.close()
        if sampAttribCsvFileInp:
            sampAttribCsvFileInp.close()
        #dump taxa summary table into CSV file
        makeFilePath(opt.predOutTaxaSummCsv)
        out = openCompressed(opt.predOutTaxaSummCsv,"w")
        flds.append("cnt")
        w = csv.DictWriter(out, fieldnames=flds, restval='null',dialect='excel-tab')
        w.writerow(dict([(fld,fld) for fld in flds]))
        for idscore,rec in sorted(idScoreSumm.items()):
            for fld_name in rec.keys():
                if fld_name.startswith("score_"):
                    rec[fld_name] /= rec["cnt"]
            w.writerow(rec)
        out.close()
        #finalize and write Krona file
        krona.finishCount()
        if not opt.isUndef("predOutStatsKronaEmbed"):
            #abspath(join(a,b)) works correctly even when b=='/c'
            krona.embedSource(
                    os.path.abspath(
                        os.path.join(
                            os.path.dirname(opt.predOutStatsHtml),
                            opt.predOutStatsKronaEmbed
                            )
                        )
                    )
        krona.write(htmlOut=opt.predOutStatsHtml,kronaUrl=opt.predOutStatsKronaUrl)
        
    def _loadPredictionsIntoSql(self,csvFile,isSumm=True):
        opt = self.opt    
        taxaLevels = self.getTaxaLevels()
        levNames = [ lev for lev in taxaLevels.getLevelNames("ascend") \
                if lev not in self.skipLevNamesInPredOut ]
        fieldsMap={"len":SqlField(type="integer")}
        if isSumm:
            fieldsMap["cnt"] = SqlField(type="integer")
            indices = {"names":("cnt","taxid","name","len","idscore","namescore")}
            tblName = "scaff_pred_summ"
        else:
            indices = {"names":("id","taxid","name","len","idscore","namescore")}
            tblName = "scaff_pred"
        for levName in levNames:
            fieldsMap["score_"+levName] = SqlField(type="real")
        makeFilePath(opt.predOutDbSqlite)
        db = DbSqlLite(dbpath=opt.predOutDbSqlite,strategy="exclusive_unsafe")
        db.createTableFromCsv(name=tblName,
                csvFile=csvFile,
                hasHeader=True,
                fieldsMap=fieldsMap,
                indices=indices)
        db.close()

    def _sqlReport(self,db,dbTable,levName,outCsv=None):
        if levName:
            fldGrp = "name_"+levName
            fldScore = "score_"+levName
            fldTaxid = "taxid_"+levName
        else:
            fldGrp = "name"
            fldScore = 0
            fldTaxid = "taxid"
        dbTableRep = dbTable+"_grp_"+levName
        sql = """\
                select %(fldGrp)s as clade,
                %(fldTaxid)s as taxid,
                sum(weight) as sum_weight,
                sum(cnt) as cnt_samp,
                sum(len) as len_samp,
                sum(len)/sum(cnt) as avg_len_samp,
                sum(%(fldScore)s*cnt)/sum(cnt) as avg_score
                from %(dbTable)s
                group by %(fldGrp)s
                order by len_samp desc
        """ % dict(fldGrp=fldGrp,fldScore=fldScore,fldTaxid=fldTaxid,dbTable=dbTable)
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
    
    def _graphicsReport(self,db,dbTableRep,levName,outBackend,fldRep="len_samp",maxClades=20):
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
        levNames = [ lev for lev in taxaLevels.getLevelNames("ascend") \
                if lev not in self.skipLevNamesInPredOut ]
        db = DbSqlLite(dbpath=opt.predOutDbSqlite,strategy="exclusive_unsafe")
        #rmrf(opt.predOutStatsDir)
        makeFilePath(opt.predOutStatsCsv)
        outCsv = openCompressed(opt.predOutStatsCsv,'w')
        import matplotlib
        matplotlib.use('AGG')
        from matplotlib.backends.backend_pdf import PdfPages
        outBackend = PdfPages(opt.predOutStatsPdf)
        sqlAsComment = True
        for levName in levNames+[""]: #"" is for lowest clade
            #TODO: add idscore level report
            dbTableRep = self._sqlReport(db=db,dbTable="scaff_pred_summ",levName=levName,outCsv=outCsv)
            self._graphicsReport(db=db,dbTableRep=dbTableRep,levName=levName,fldRep="sum_weight",
                    outBackend=outBackend,maxClades=20)
        outCsv.close()
        outBackend.close()
        db.close()

    def procBenchScores(self,**kw):
        """Process benchmark IMM scores and generate performance metrics.
        Parameters are taken from self.opt.
        @param dbBench path of ImmClassifierBenchmark instance that was used
        in scoring. Metadata from this instance will be used here.
        @param outScoreComb File with ImmScores object
        @param outIdScoreMeta Metadata map for each score ID in outScoreComb
        @param benchOutDir Output directory for benchmarking results
        @param benchOutDbSqlite Output file with metrics in SQLite format
        """
        opt = self.opt
        makedir(opt.benchOutDir)
        makeFilePath(opt.benchOutDbSqlite)
        #set i_lev_excl to this value when no exclusion was done
        iLevNoExc = 0
        db = DbSqlLite(dbpath=opt.benchOutDbSqlite,strategy="exclusive_unsafe")
        debugShortcut = False
        if not debugShortcut:
            storesDb = self.prepBenchSql(db=db)
            tblLevels = storesDb.storeDbLev.tblLevels
            tblNames = storesDb.storeDbNode.tblNames
            tblModNames = storesDb.tblModNames
        else:
            tblLevels = "txlv"
            tblNames = "taxa_names"
            tblModNames = "model_names"
        self.procBenchScoreMatr(db=db,iLevNoExc=iLevNoExc)
        sqlBenchMetr = ImmClassifierBenchMetricsSql(db=db,
                taxaLevelsTbl=tblLevels,
                taxaNamesTbl=tblNames,
                modelNamesTbl=tblModNames,
                iLevNoExc=iLevNoExc)
        sqlBenchMetr.makeMetrics(nLevTestModelsMin=opt.benchNLevTestModelsMin,
                outDir=opt.benchOutDir,
                comment=("fragment length %s" % (opt.dbBenchFragLen,)) \
                        if opt.dbBenchFragLen \
                        else None)

    def prepBenchSql(self,db):
        """Prepare lookup tables in benchmark SQL database.
        @param db Sql object"""
        opt = self.opt
        idScoreMeta = loadObj(opt.outIdScoreMeta)
        tblModNames="model_names"
        db.createTableFromKeyVal(tblModNames,
                records=(
                    (key,val["name"]) for (key,val) in idScoreMeta.items()
                    ),
                keyName="id",
                keyType="char({})".format(maxIdLen),
                valName="name",
                valType="text"
                )
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        print "DEBUG: Loaded taxa tree and levels"
        tableSfx = ""
        storeDbNode = NodeStorageDb(db=db,tableSfx=tableSfx)
        storeDbNode.saveName(taxaTree)
        storeDbLev = LevelsStorageDb(db=db,tableSfx=tableSfx)
        storeDbLev.save(taxaLevels)
        return Struct(storeDbLev=storeDbLev,
                storeDbNode=storeDbNode,
                tblModNames=tblModNames)

    def procBenchScoreMatr(self,db,iLevNoExc,**kw):
        """Process benchmark IMM scores and generate (test,pred) pairs in SQL.
        @param db Sql object

        Other parameters are taken from self.opt.
        @param dbBench path of ImmClassifierBenchmark instance that was used
        in scoring. Metadata from this instance will be used here.
        @param outScoreComb File with ImmScores object
        @param outIdScoreMeta Metadata map for each score ID in outScoreComb
        @param benchOutCsv Output file with metrics in CSV format
        @param benchOutDbSqlite Output file with metrics in SQLite format
        @param benchProcSubtrees Filter samples by taxonomy IDs in this file
        """
        opt = self.opt
        immScores = openImmScores(opt,fileName=opt.outScoreComb,mode="r")
        sc = immScores.getData()
        idScoreMeta = loadObj(opt.outIdScoreMeta)
        #assume idScore are idImm:
        idImm = sc.idScore[:]
        kind = immScores.getKind()
        assert kind == "ImmScoresDenseMatrix","Matrix Score dataset is required for benchmark"
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        print "DEBUG: Loaded taxa tree and levels"
        eukRoot = taxaTree.getNode(eukTaxid)
        levNames = taxaLevels.getLevelNames("ascend")
        levIdNames = taxaLevels.getLevelIdForNames(levNames)
        levPosSuperking = levNames.index("superkingdom")
        levPosSpecies = levNames.index("species")
        levPosGenus = levNames.index("genus")
        max_linn_levid = taxaLevels.getLevelId("order")
        min_linn_levid = taxaLevels.getLinnLevelIdRange()[0]
        #this value as prediction result should be treated as
        #always incorrect when computing performance metrics
        wrongTaxid = self.rejTaxid - 1
        sampTaxidsMap = None
        if opt.benchTaxidsMap is not None:
            raise ValueError("Depricated option --bench-taxids-map. "+\
                    "Define --bench-id-seq-db-to-id-score-remapping instead")
        if opt.benchIdSeqDbToIdScoreRemapping is not None:
            idSeqDbToIdScoreRemapping = \
                    load_config_json(opt.benchIdSeqDbToIdScoreRemapping)
        else:
            idSeqDbToIdScoreRemapping = None
        benchDb = ImmClassifierBenchmark.open(opt.dbBench,"r") 
        idSampToIdScore = benchDb.mapIdFragToIdScore(idScoreMeta=idScoreMeta,
                idSeqDbToIdScoreRemapping=idSeqDbToIdScoreRemapping)
        #sampTaxids = set()
        #for iS,idS in enumerate(sc.idSamp[:]):
        #    taxid,itS = idS.strip().split('_')
        #    taxid,itS = int(taxid),int(itS)
        #    sampTaxids.add(taxid)

        def _mapModelIds(taxaTree,idImm,idScoreMeta):
            """Map model ids to the taxonomy tree and store for each node the
            number of immediate children that have models in their subtrees.
            This will be used to include into the performance counters only
            those samples for which a correct model exists at a given exclusion
            level"""
            root = taxaTree.getRootNode()
            root.setAttribute("tmpHasMod",0)
            root.setAttribute("tmpHasLeafMod",0)
            for id_imm in idImm:
                meta = idScoreMeta[id_imm]
                node = taxaTree.getNode(meta["taxid"])
                node.tmpHasMod = 1
                node.tmpHasLeafMod = int(meta["is_leaf"]==True)
            root.setTotal("tmpHasMod","tmpCntMod")
            #root.setReduction(extractor=lambda n: 0,
            #        dstAttr="tmpChWMod",
            #        childExtractor=lambda n: n.tmpCntMod > 0)
            root.setTotal("tmpHasLeafMod","tmpCntLeafMod")

        _mapModelIds(taxaTree=taxaTree,idImm=idImm,idScoreMeta=idScoreMeta)
        print "DEBUG: mapped models to tree"
        
        def _mapBenchProcSubtreesToTree(taxaTree,benchProcSubtreesFiles):
            """Mark those nodes for which to process samples"""
            if benchProcSubtreesFiles:
                taxids = set()
                for benchProcSubtreesFile in benchProcSubtreesFiles:
                    with openCompressed(benchProcSubtreesFile,"r") as inp:
                        for line in inp:
                            taxids.add(int(line.strip().split('\t')[0]))
                for taxid in taxids:
                    try:
                        node = taxaTree.getNode(taxid)
                        node.setAttribute("tmpProcSamp",1)
                    except KeyError:
                        print "DEBUG: benchProcSubtree taxid %s not found in a tree" % (taxid,)
        _mapBenchProcSubtreesToTree(taxaTree=taxaTree,
                benchProcSubtreesFiles=opt.benchProcSubtrees)
        print "DEBUG: mapped --bench-proc-subtrees to a tree"

        db.ddl("""
        create table bench_samp 
        ( 
        i_samp integer NOT NULL,
        i_lev_exc integer NOT NULL,
        i_lev_per integer NOT NULL,
        id_mod_test character({id_mod_len}),
        id_mod_pred character({id_mod_len}),
        taxid_lev_test_bot integer,
        taxid_lev_pred_bot integer,
        taxid_lev_test integer,
        taxid_lev_pred integer,
        taxid_superking integer,
        n_lev_test_models integer,
        score real
        )
        """.format(id_mod_len=maxIdLen),
        dropList=["table bench_samp"])
        inserter = db.makeBulkInserter(table="bench_samp")
        def _outTestRec(iS,
                levIdFS,
                levIdFP,
                idModS,
                idModP,
                nodS,
                nodP,
                nodLinFS_FP,
                nodLinFS,
                nodLinFP,
                nodSuperkingS,
                scoP,
                withExcl,
                n_lev_test_models,
                inserter=inserter,levIdNames=levIdNames,
                wrongTaxid=wrongTaxid,
                iLevNoExc=iLevNoExc):
            rec = (iS,
                    levIdFS,
                    levIdFP,
                    idModS,
                    idModP,
                    nodS.id,
                    nodP.id,
                    nodLinFS_FP.id,
                    nodLinFP.id if nodLinFP is not None else wrongTaxid,
                    nodSuperkingS.id,
                    n_lev_test_models,
                    scoP)
            inserter(rec)
        #If we select columns directly from hdf file with row-optimized
        #chunkshape, the performance is horrible and we spend 100% CPU in
        #system mode calls. For now, we load entire array into numpy.
        #@todo add copy to column-optimized chunkshape to make the final 
        #hdf5 (per pytables' author recipe), then select from file directly.
        score = sc.score[:]
        print "DEBUG: Loaded score matrix of total elements %i" % (score.shape[0]*score.shape[1],)
        #rows are models, columns are samples
        #samples ID is UUID; samples were organized into SampDbFasta indexed by idSeqDb
        #that is in turn listed in idScoreMeta["idScore"]["seq_db_id"]
        #we cache variables that depend only on true sample taxid or predicted sample
        #taxid because the former go in series, and the later often too (when predictions
        #are correct, for instance)
        midS = None #model id sample (true value)
        tidS = None #taxid sample (true value)
        nodS = None #node(taxid sample)
        linS = None #lineage sample
        linSetS = None #lineage sample as set
        linFS = None #lineage fixed list for sample
        midP = None #model id predicted
        tidP = None #taxid predicted
        nodP = None #node(taxid predicted)
        linP = None #lineage predicted
        fill = None #None #"up-down" "up"
        assert fill is None,"The consequences of using fill for fixed lineage are not fully understood yet"
        cntRej = defdict(int)
        cntPass = defdict(int)
        iSel = sorted(sampleRatioWOR(n.arange(len(sc.idSamp)).tolist(),0.1))
        #for iS,idS in it.izip(iSel,sc.idSamp[iSel]):
        for iS,idS in enumerate(sc.idSamp):
            idMod = idSampToIdScore[idS]
            if midS is None or midS != idMod:
                midS = idMod
                modMeta = idScoreMeta[midS]
                tidS = modMeta["taxid"]
                nodS = taxaTree.getNode(tidS)
                linS = nodS.lineage()
                linSetS = set(linS)
                linFS = taxaLevels.lineageFixedList(nodS,null=None,format="node",fill=fill)
            #model can be attached to any level in the taxonomic hierarchy
            if nodS.tmpHasLeafMod and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #if nodS.tmpHasLeafMod and linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #if linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #if nodS.tmpHasLeafMod and taxidOrig not in [83129,484896,693582,490912,490913] and linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #lytic
            #if taxidOrig in [83129,484896,693582,490912,490913] and linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #lysogenic
            #if taxidOrig in [151599,397353,373126,186152,198539] and linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #if taxidOrig in [151536,10730,373126,186152,198539] and linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #if not nodS.tmpHasLeafMod and linFS[levPosSpecies] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
            #if linFS[levPosGenus] and (not opt.benchProcSubtrees or hasattr(nodS,"tmpProcSamp")):
                scoreP = score[:,iS]
                ordP = (-(scoreP)).argsort()
                startOrdP = 0
                nodSuperkingS = linFS[levPosSuperking]
                #nodLinFS_Prev = None
                #iLevFS is exclusion level
                #we additionally prepend the node itself to process w/o exclusion
                levTasks = [ (True,iLevFS,nodLinFS) for iLevFS,nodLinFS in enumerate(linFS[:-1]) ]
                levTasks = [ (False,iLevNoExc,nodS) ] + levTasks
                    
                #for iLevFS,nodLinFS in enumerate(linFS[:-1]):
                #    if nodLinFS is not None:
                #        levTasks = [ (False,iLevFS,nodLinFS) ] + levTasks
                #        break
                
                for withExcl,iLevFS,nodLinFS in levTasks:
                    outputMetrForLineage = False
                    if nodLinFS is not None: 
                        #The condition above checks that the exclusion level is present in lineage 
                        #of a true node. Clearly, we cannot speak of excluding related models 
                        #at a certain level if the test sample does not have that level defined.
                        #So, the above omits such records from consideration regardless of the
                        #predicted taxa.
                        for iOrdP,vOrdP in enumerate(ordP,startOrdP):
                            scoP = float(scoreP[vOrdP])
                            midP = idImm[vOrdP]
                            tidP = idScoreMeta[midP]["taxid"]
                            nodP = taxaTree.getNode(tidP)
                            #TMP:
                            if False and not taxaLevels.isNodeInLinnLevelRange(nodP,\
                                    min_linn_levid,max_linn_levid):
                                break
                            if nodP.tmpHasLeafMod: #excludes mixed models
                                if withExcl:
                                    if \
                                        ((nodP not in linSetS \
                                        and not nodP.isUnder(nodLinFS))):
                                            outputMetrForLineage = True
                                            break
                                else:
                                    outputMetrForLineage = True
                                    break

                        if outputMetrForLineage:
                            cntPass[(withExcl,levIdNames[iLevFS])] += 1
                        else:
                            cntRej[(withExcl,levIdNames[iLevFS])] += 1

                        if outputMetrForLineage:
                            linFP = taxaLevels.lineageFixedList(nodP,null=None,format="node",fill=fill)
                            #climb the levels and output (test,pred) pairs
                            #iLevFP is prediction level > iLevFS or == iLevFS if not withExcl.
                            if withExcl:
                                iLevFP_start = iLevFS + 1
                            else:
                                iLevFP_start = iLevFS
                            #Check that sample is not the only child with model at the lowest
                            #define prediction level (or that there is a least one model if we are not excluding clades)
                            #Phymm Methods for unexplained reason says they required 
                            #"two others", which must mean 3 total, so we also save the 
                            #number as part of the record to filter on it later when 
                            #computing the metrics.
                            #Note that checking instead for sisters at each prediction level in the lineage
                            #would be wrong because it would be equal to raising the exclusion level for some
                            #nodes leading to cases when a family prediction level has lower average metrics
                            #than a genus level.
                            n_lev_test_models = 0
                            if withExcl:
                                for iLevFP,nodLinFP in enumerate(linFP[iLevFP_start:],iLevFP_start):
                                    nodLinFS_FP = linFS[iLevFP] #true node at prediction level
                                    if nodLinFS_FP is not None: 
                                        n_lev_test_models = nodLinFS_FP.tmpCntLeafMod - nodLinFS.tmpCntLeafMod
                                        break
                            else:
                                n_lev_test_models = nodLinFS.tmpCntLeafMod
                            if n_lev_test_models >= 1:
                                for iLevFP,nodLinFP in enumerate(linFP[iLevFP_start:],iLevFP_start):
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

                                        #if nodLinFS_FP.id == 1081:
                                        #    pdb.set_trace()
                                        _outTestRec(iS,
                                                levIdNames[iLevFS] if withExcl else iLevNoExc,
                                                levIdNames[iLevFP],
                                                midS,
                                                midP,
                                                nodS,
                                                nodP,
                                                nodLinFS_FP,
                                                nodLinFS,
                                                nodLinFP,
                                                nodSuperkingS,
                                                scoP,
                                                withExcl,
                                                n_lev_test_models)
                                    
                                if not withExcl:
                                    #insert also the record for the test with the predicted
                                    #nodes itself when no exclusion is done
                                    _outTestRec(iS,
                                            iLevNoExc,
                                            iLevNoExc,
                                            midS,
                                            midP,
                                            nodS,
                                            nodP,
                                            nodS,
                                            nodS,
                                            nodP,
                                            nodSuperkingS,
                                            scoP,
                                            withExcl,
                                            n_lev_test_models)

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
        "i_samp",
        "i_lev_exc",
        "i_lev_per",
        "id_mod_test",
        "id_mod_pred",
        "taxid_lev_test",
        "taxid_lev_pred",
        "taxid_lev_test_bot",
        "taxid_lev_pred_bot",
        "taxid_superking",
        "n_lev_test_models"
        ])
        cntPass = dict(cntPass)
        cntRej = dict(cntRej)
        print "cntRej=%s\tcntPass=%s" % (sorted(cntRej.items()),sorted(cntPass.items()))

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmClassifierApp)

