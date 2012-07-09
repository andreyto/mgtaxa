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


class MGTOptions(Options):
    def __init__(self,*l,**kw):
        Options.__init__(self,*l,**kw)
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
        self.toolGpgKeyName = "mgtaxa@jcvi.org"
        self.tmpDir = os.environ["MGT_TMP"]
        self.dataDir = os.environ["MGT_DATA"]
        self.ncbiDbDir = os.environ["MGT_NCBI_DB"]
        self.srcDir = "/home/"+os.environ["USER"]+"/work/mgtaxa"
        def _pdata(dir):
            return pjoin(self.dataDir,dir)
        def _pdata_ncbi(dir):
            return pjoin(self.ncbiDbDir,dir)
        self.homeDir = os.environ["MGT_HOME"]
        self.instDir = self.homeDir
        self.binDir = pjoin(self.instDir,"bin")
        self.confDir = pjoin(self.instDir,"etc")
        self.envRc = os.environ["MGT_RC"]
        self.testDataDir = os.path.join(self.srcDir,"test_data")
        self.testRunDir = os.path.join(self.srcDir,"test_run")
        self.refSeqDataDir = _pdata_ncbi("refseq")
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
        self._setTaxaFileNames(_pdata_ncbi('taxonomy'),sfx="")
        self._setTaxaFileNames(_pdata_ncbi('taxonomy.new'),sfx="New")
        self._setTaxaFileNames(pjoin(self.testDataDir,"taxonomy"),sfx="Test")
        self.taxaTreeTablePrefix = "txtr"
        self.taxaTreeTableSfxMain = ""
        self.taxaLevelsTablePrefix = "txlv"
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
                lrmSiteOptions=r'-l memory=2000M -l arch="lx*64"',
                lrmUserOptions='-P 0413',
                maxQueued = 1500,
                envRc = self.envRc,
                #[qsub|makeflow]
                batchBackend="qsub")
        self.batchRunTerminator = deepcopy(self.batchRun)
        self.batchRunTerminator.lrmSiteOptions=r'-l memory=200M -l fast -l arch="lx*64"'
        #terminator only makes sense outside of makeflow
        self.batchRunTerminator.batchBackend="qsub"
        self.app = Options(
                #override App.opt.runMode value if not "default" here. Other choices are "inproc"
                runMode = "default",
                #extra default arguments to python that executes the App script
                extraPyArgs = None 
                )
        self.genomeTools = Options(
                exe=pjoin(os.environ["INSTMACH"],"bin/gt"),
                sketchConf=pjoin(self.confDir,"gt.sketch.default.style")
                )
        self.glimmer3 = Options(
                topDir=os.environ["INSTMACH"]
                )
        self.glimmer3.immBuildExe=pjoin(self.glimmer3.topDir,"bin","mgt-glm-build-icm")
        self.glimmer3.immScoreExe=pjoin(self.glimmer3.topDir,"bin","mgt-glm-simple-score")
        self.icm = Options(
                icmDb = pjoin(self.dataDir,"icm-refseq"),
                predScoreRescaleModel = {
                    "lreg" : {
                        "coef" : {
                            "(Intercept)":20.9733729119449,
                            "score_b":15.3156079336066,
                            "len_samp.pre":-77.8559030330633,
                            "i_lev_per.pre":0.0382037195511002,
                            "score_b:len_samp.pre":-43.8872972190474
                            }
                        }
                    }
                )

        self.krona = Options(
                topDir=pjoin(os.environ["INSTMACH"],"krona")
                )
        self.krona.cmdBase=["perl","-I",pjoin(self.krona.topDir,"lib")]
        self.krona.xmlToChart=self.krona.cmdBase+[pjoin(self.krona.topDir,"scripts/ImportXML.pl")]


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

