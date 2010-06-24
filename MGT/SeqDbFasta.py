"""A simple "database" of multi-FASTA files."""

from MGT.DirStore import *
from MGT.Taxa import *
from MGT.FastaIO import *
from MGT.FeatCommon import *

class SeqDbFasta(DirStore):
    """An interface to a collection of FASTA files.
    Each file is named with some ID, such as taxonomy ID.
    Methods for creating the DB grouped by taxonomy and streaming
    based on lists of IDs are provided.
    """

    fastaSfx = ".fasta.gz"

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def setTaxaTree(self,taxaTree):
        self.taxaTree = taxaTree

    def importByTaxa(self,inpFiles):
        taxaTree = self.getTaxaTree()
        splitFastaFilesByTaxa(inSeqs=inpFiles,taxaTree=taxaTree,giToTaxa=None,outDir=self.getPath())

    def loadTaxaList(self):
        self.taxaList = n.asarray([ int(f) for f in self.fileNames("*"+self.fastaSfx,sfxStrip=self.fastaSfx) ])

    def getTaxaList(self):
        if not hasattr(self,"taxaList"):
            self.loadTaxaList()
        return self.taxaList

    getIdList = getTaxaList

    def getFilePathById(self,id):
        return self.getFilePath("%s%s" % (id,self.fastaSfx))

    def fastaReader(self,id):
        return FastaReader(self.getFilePathById(id))

    def writeFasta(self,ids,out):
        for id in ids:
            reader = self.fastaReader(id)
            for rec in reader.records():
                out.write(rec.header())
                for line in rec.seqLines():
                    out.write(line)
            reader.close()

    def writeFastaBothStrains(self,ids,out):
        fsout = None # create after lineLen is known
        for id in ids:
            reader = self.fastaReader(id)
            for rec in reader.records():
                header = rec.header()
                seq = rec.sequence()
                if len(seq) > 0:
                    if fsout is None: 
                        fsout = FastaWriter(out,lineLen=reader.lineLen())
                    fsout.record(header=header,sequence=seq)
                    fsout.record(header=header,sequence=revCompl(seq))
            reader.close()

