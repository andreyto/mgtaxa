from MGT.Util import Options
import os

class MGTOptions:
    def __init__(self):
        self.debug = 1
        # Use modified 'fastacmd' binary that replaces nucleotide
        # content with repeated 'x<gi>' string. This should be set
        # before TaxaCollector.rebuild() is called. It will activate
        # several assertions throught the code up to the k-mer generation
        # stage, with the goal of making sure that we do not mix up sequence
        # content from different taxonomy IDs when picking trainig/testing
        # samples etc.
        self.debugFakeSequence = False
        self.tmpDir = "/usr/local/scratch/atovtchi"
        self.dataDir = "/home/atovtchi/work/mgtdata"
        self.taxaPickled = 'taxa.pkl.gz'
        # Max length to store for a full header. It is stored in a separate table,
        # so size does not matter that much
        self.fastaHdrSqlLen = 2048
        self.blastDataDir = 'blast'
        self.collectTaxaGiFile = os.path.abspath('mgt_coll_taxa.gi')
        self.seqGiIdFile = os.path.abspath('mgt_seq.giid')
        #self.selDumpFile = 'mgt_sel.csv'
        self.selFastaFile = 'mgt_sel.fasta.gz'
        self.blastSelAlias = 'mgt'
        self.taxaDataDir = os.path.join(self.dataDir,'taxonomy')
        self.taxaCatFile = os.path.join(self.taxaDataDir,'categories.dmp')
        self.taxaGiFile = os.path.join(self.taxaDataDir,'gi_taxid_nucl.dmp.gz') 
        self.taxaDumpDir = self.taxaDataDir
        self.taxaNodesFile = os.path.join(self.taxaDumpDir,'nodes.dmp')
        self.taxaDivisionFile = os.path.join(self.taxaDumpDir,'division.dmp')
        self.taxaNamesFile = os.path.join(self.taxaDumpDir,'names.dmp')
        self.taxaTreeTablePrefix = "txtr"
        self.taxaTreeTableSfxMain = ""
        #self.kmerTxtDir = os.environ['PHYLA_KMERS']
        #self.kmerTestFile = os.path.join(self.kmerTxtDir,'6mers_1K.gz')
        self.hdfSeqFile = 'seq.hdf'
        self.hdfSeqGroup = '/seq'
        #This will only hold index that references sequence from hdfSeqFile
        self.hdfActSeqFile = 'act_seq.hdf'
        self.hdfActSeqInd = '/ind'
        #This will hold sample index that references both active sequence index
        #and sequence data. We create a separate HDF file for each sample chunkLen
        self.hdfSampGroup = '/samp'

        sampSel = Options(
                all=Options(
                    longSeq = 100000,
                    prod=Options(
                        train=Options(min=40,max=82000) #40 35000
                        ),
                    test=Options(
                        test=Options(min=5,ratio=0.3) #5
                        )
                    )
                )
        test_train = sampSel.all.prod.train.copy()
        # this will result in a fewer number of trainable nodes
        # because selecting at least test.min samples usually
        # subtracts more than that from available training samples
        # (since we select for testing entire genus or species nodes)
        test_train.min -= sampSel.all.test.test.min
        sampSel.all.test.train = test_train
        vir = sampSel.all.copy()
        vir.longSeq = 0
        vir.prod.train.min = 30 #30
        vir.test.train.min = vir.prod.train.min - vir.test.test.min
        sampSel.vir = vir
        sampSel.freeze()
        self.sampSel = sampSel
        self.maxTestSampLen = 1000
        self.minTestSampLen = 1000
        self.labelTopNodeId = 10239
        self.labelLevel = "family"
        self.sampLen = 1000
        self.kmerLen = 8
        self.maxTestSampPerTaxa = 30
        #self.kmerRepr = "Frequences"
        #self.kmerRepr = "Bits"
        self.kmerRepr = "Sequences"
        self.sampNamePrefix = "samp_1k"
        self.predictorDir = os.path.join(self.tmpDir,"pred-1k-vfam2-1")
        self.hdfTestFile = 'test.hdf'
        self.svmTestFile = 'test.svm'
        self.svmTrainFile = 'train.svm'
        self.predictorTable = "pred"
        self.batchRun = Options(
                PROJECT_CODE = 600005,
                MEM = 1000,
                ARCH = "lx26-eon64",
                maxQueued = 50)

options = MGTOptions()

