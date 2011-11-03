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

    def importByTaxa(self,inpFiles,filt=None):
        taxaTree = self.getTaxaTree()
        splitFastaFilesByTaxa(inSeqs=inpFiles,taxaTree=taxaTree,giToTaxa=None,outDir=self.getPath(),filt=filt)

    def loadTaxaList(self):
        self.taxaList = n.asarray([ int(f) for f in self.fileNames("*"+self.fastaSfx,sfxStrip=self.fastaSfx) ])

    def getTaxaList(self):
        if not hasattr(self,"taxaList"):
            self.loadTaxaList()
        return self.taxaList

    getIdList = getTaxaList

    def getFileBaseById(self,id):
        return "%s%s" % (id,self.fastaSfx)

    def getFilePathById(self,id):
        return self.getFilePath(self.getFileBaseById(id))

    def updateMetaDataById(self,id):
        meta = Struct()
        meta["seqLengths"] = self.computeSeqLengths(id)
        self.saveFileMetaData(meta,self.getFileBaseById(id))

    def finById(self,id):
        """Finalize creation of a new ID in SeqDb after using methods such as importByTaxa.
        This creates meta-data and possibly does other things.
        @todo make it also compress the main data object if necessary"""
        self.updateMetaDataById(id)

    def delById(self,id):
        os.remove(self.getFilePathById(id))

    def fastaReader(self,id):
        """Return FastaReader to read from the DB for a given ID"""
        return FastaReader(self.getFilePathById(id))

    def computeSeqLengths(self,id):
        """Compute a numpy recarray with fields ("id","len") for a given DB ID
        directly from sequence data.
        Currently has to read all sequences for id - works through FastaReader"""
        return fastaLengths(self.getFilePathById(id))

    def seqLengths(self,id):
        """Return a numpy recarray with fields ("id","len") for a given DB ID.
        Loads from meta data that has to be pre-computed."""
        meta = self.loadFileMetaData(self.getFileBaseById(id))
        return meta["seqLengths"]
    
    def fastaWriter(self,id,lineLen=None,mode="w"):
        """Return FastaWriter to write INTO the DB for a given ID.
        @param id ID of the record
        @param lineLen Length of FASTA lines
        @param mode String with the same semantics as 'mode' supplied to built-in open() method
        """
        return FastaWriter(out=self.getFilePathById(id),lineLen=lineLen,mode=mode)

    def writeFasta(self,ids,out):
        """Write FASTA stream for given sequence of IDs from DB into output file-like object"""
        for id in ids:
            reader = self.fastaReader(id)
            for rec in reader.records():
                out.write(rec.header())
                for line in rec.seqLines():
                    out.write(line)
            reader.close()

    def balanceLengths(self,ids,maxLen):
        totalLength = sum([ self.seqLengths(id)["len"].sum() for id in ids ])
        lenRatio = float(maxLen)/totalLength
        return lenRatio

    def writeFastaBothStrands(self,ids,out,maxLen=None):
        totalWritten = 0
        fsout = None # create after lineLen is known
        if maxLen is not None:
            #the actual amount of sequence written will be double that
            #because of rev-compl
            lenRatio = self.balanceLengths(ids,maxLen)
            if lenRatio >= maxLen:
                lenRatio = None
        else:
            lenRatio = None
        for id in ids:
            reader = self.fastaReader(id)
            for rec in reader.records():
                header = rec.header()
                seq = rec.sequence()
                seqLen = len(seq)
                if lenRatio is not None:
                    #@todo: subsampling each sequence is a questionable
                    #strategy when we are dealing with already short 
                    #sequences from a WGS project.
                    #at the opposing end, long sequences should be 
                    #subsampled as a series of random fragments each
                    #
                    #This gets a substring of a given length starting at a 
                    #random position within the full sequence
                    seq = SubSamplerRandomStart(int(len(seq)*lenRatio))(seq)
                if len(seq) > 0:
                    if fsout is None: 
                        fsout = FastaWriter(out,lineLen=reader.lineLen())
                    fsout.record(header=header,sequence=seq)
                    fsout.record(header=header,sequence=revCompl(seq))
                    print "DEBUG:     id=%s seqLen=%s seqLenWritten=%s" % (id,seqLen,len(seq))
                    totalWritten += len(seq)
            reader.close()
        print "DEBUG: id=%s maxLen=%s lenRatio=%s totalWritten=%s" % (id,maxLen,lenRatio,totalWritten)

