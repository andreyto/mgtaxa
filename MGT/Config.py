### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Options import *
from MGT.TaxaConst import *
import os
pjoin = os.path.join
pabs = os.path.abspath


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
        self.toolName = "MGT"
        self.toolEmail = "atovtchi@jcvi.org"
        self.tmpDir = "/usr/local/scratch/atovtchi"
        self.dataDir = "/home/atovtchi/work/mgtdata"
        self.srcDir = "/home/atovtchi/work/mgtaxa"
        def _pdata(dir):
            return pjoin(self.dataDir,dir)
        self.homeDir = os.environ["MGT_HOME"]
        self.instDir = self.homeDir
        self.binDir = pjoin(self.instDir,"bin")
        self.confDir = pjoin(self.instDir,"etc")
        self.envRc = os.environ["MGT_RC"]
        self.testDataDir = os.path.join(self.srcDir,"test_data")
        self.refSeqDataDir = _pdata("refseq")
        # Max length to store for a full header. It is stored in a separate table,
        # so size does not matter that much
        self.fastaHdrSqlLen = 2048
        ## Drop sequence records with these NCBI divid's at the earliest stage of DB construction
        ## Currently this is unknown stuff + large creatures (that still leaves plants because
        ## they include algae)
        self.dividDrop= (7,8,11) + (2,5,6,10)
        ## Drop sequence records in subtrees (top included) under these NCBI taxid's at the earliest stage of DB construction
        self.taxidDrop = (embryophytaTaxid,chordataTaxid)
        self.blastDataDir = 'blast'
        self.collectTaxaGiFile = _pdata('mgt_coll_taxa.gi')
        self.seqGiIdFile = _pdata('mgt_seq.giid')
        #self.selDumpFile = 'mgt_sel.csv'
        self.selFastaFile = 'mgt_sel.fasta.gz'
        self.blastSelAlias = 'mgt'
        self._setTaxaFileNames(_pdata('taxonomy'),sfx="")
        self._setTaxaFileNames(_pdata('taxonomy.new'),sfx="New")
        self.taxaTreeTablePrefix = "txtr"
        self.taxaTreeTableSfxMain = ""
        #self.kmerTxtDir = os.environ['PHYLA_KMERS']
        #self.kmerTestFile = os.path.join(self.kmerTxtDir,'6mers_1K.gz')
        self.hdfSeqFile = _pdata('seq.hdf')
        self.hdfSeqGroup = '/seq'
        #This will only hold index that references sequence from hdfSeqFile
        self.hdfActSeqFile = _pdata('act_seq.hdf')
        self.hdfActSeqInd = '/ind'
        #This will hold sample index that references both active sequence index
        #and sequence data. We create a separate HDF file for each sample chunkLen
        self.hdfSampGroup = '/samp'

        sampSel = Options(
                all=Options(
                    longSeq = 100000,
                    prod=Options(
                        train=Options(min=40,max=35000) #40 35000
                        ),
                    test=Options(
                        test=Options(min=5,ratio=0.1) #5 0.3
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
        #vir.prod.train.min = 30 #30
        #vir.test.train.min = vir.prod.train.min - vir.test.test.min
        sampSel.vir = vir
        sampSel.freeze()
        self.sampSel = sampSel
        self.maxTestSampLen = 5000
        self.minTestSampLen = 5000
        self.labelTopNodeId = 1 #35237 #10239 - all vir, 35237 - dsDNA vir
        #self.labelLevel = "superkingdom"
        #self.labelLevel = "family"
        #self.labelLevel = "order"
        self.sampLen = 5000
        self.kmerLen = 8
        self.maxTestSampPerTaxa = 100
        #self.kmerRepr = "Frequences"
        self.kmerRepr = "Bits"
        #self.kmerRepr = "Sequences"
        self.sampNamePrefix = "samp_5000"
        self.hdfTestFile = 'test.hdf'
        self.predictorTable = "pred"
        self.batchRun = Options(
                PROJECT_CODE = 600005,
                MEM = 2000,
                ARCH = "lx*64",
                maxQueued = 50,
                LENGTH = "medium", #fast
                ENVRC = self.envRc)
        self.app = Options(
                runMode = "default"#override App.opt.runMode value if not "default" here. Other choices are "inproc"
                )
        self.genomeTools = Options(
                exe=pjoin(os.environ["INSTMACH"],"app/gt/gt"),
                sketchConf=pjoin(self.confDir,"gt.sketch.default.style")
                )

    def _setTaxaFileNames(self,topDir,sfx=""):
        setattr(self,'taxaDataDir'+sfx,topDir)
        taxaDataDir = getattr(self,'taxaDataDir'+sfx) 
        setattr(self,'taxaPickled'+sfx,os.path.join(taxaDataDir,'gi_taxid.pkl.gz'))
        setattr(self,'taxaCatFile'+sfx,os.path.join(taxaDataDir,'categories.dmp'))
        setattr(self,'taxaGiFiles'+sfx,[ os.path.join(taxaDataDir,'gi_taxid_%s.dmp.gz' % moltype) 
                for moltype in ("nucl","prot") ])
        setattr(self,'taxaDumpDir'+sfx,taxaDataDir)
        taxaDumpDir = getattr(self,'taxaDumpDir'+sfx) 
        setattr(self,'taxaNodesFile'+sfx,os.path.join(taxaDumpDir,'nodes.dmp'))
        setattr(self,'taxaDivisionFile'+sfx,os.path.join(taxaDumpDir,'division.dmp'))
        setattr(self,'taxaNamesFile'+sfx,os.path.join(taxaDumpDir,'names.dmp'))

    def getPyCmd(self,cmd):
        return "python %s.py" % pjoin(self.binDir,cmd)

options = MGTOptions()

