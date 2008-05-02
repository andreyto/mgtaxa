import os

class Options:
    def __init__(self):
        self.debug = 1
        self.tmpDir = "/usr/local/scratch/atovtchi"
        self.taxaPickled = 'taxa.pkl.gz'
        # Max length to store for a full header. It is stored in a separate table,
        # so size does not matter that much
        self.fastaHdrSqlLen = 2048
        self.blastDataDir = 'blast'
        self.collectTaxaGiFile = os.path.abspath('mgt_coll_taxa.gi')
        #self.selDumpFile = 'mgt_sel.csv'
        self.selFastaFile = 'mgt_sel.fasta.gz'
        self.blastSelAlias = 'mgt'
        self.taxaDataDir = 'taxonomy'
        self.taxaCatFile = os.path.join(self.taxaDataDir,'categories.dmp')
        self.taxaGiFile = os.path.join(self.taxaDataDir,'gi_taxid_nucl.dmp.gz') 
        self.taxaDumpDir = self.taxaDataDir
        self.taxaNodesFile = os.path.join(self.taxaDumpDir,'nodes.dmp')
        self.taxaDivisionFile = os.path.join(self.taxaDumpDir,'division.dmp')
        self.taxaNamesFile = os.path.join(self.taxaDumpDir,'names.dmp')
        self.taxaTreeTablePrefix = "txtr"
        #self.kmerTxtDir = os.environ['PHYLA_KMERS']
        #self.kmerTestFile = os.path.join(self.kmerTxtDir,'6mers_1K.gz')
        self.hdfCollTaxaFile = 'collTaxa.hdf'
        self.hdfCollTaxaGroup = 'seq'


options = Options()
