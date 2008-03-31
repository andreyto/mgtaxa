import os

class Options:
    def __init__(self):
        self.taxaPickled = 'taxa.pkl.gz'
        self.fastaHdrSqlLen = 40
        self.blastDataDir = 'blast'
        self.selGiFile = 'phyla_sel.gi'
        self.selDumpFile = 'phyla_sel.csv'
        self.selFastaFile = 'phyla_sel.fasta.gz'
        self.srcDbNameAlias = 'phyla'
        self.taxaDataDir = 'taxonomy'
        self.taxaCatFile = os.path.join(self.taxaDataDir,'categories.dmp')
        self.taxaGiFile = os.path.join(self.taxaDataDir,'gi_taxid_nucl.dmp.gz') 
        self.taxaDumpDir = self.taxaDataDir
        self.taxaNodesFile = os.path.join(self.taxaDumpDir,'nodes.dmp')
        self.taxaDivisionFile = os.path.join(self.taxaDumpDir,'division.dmp')
        self.taxaNamesFile = os.path.join(self.taxaDumpDir,'names.dmp')
        #self.kmerTxtDir = os.environ['PHYLA_KMERS']
        #self.kmerTestFile = os.path.join(self.kmerTxtDir,'6mers_1K.gz')


options = Options()
