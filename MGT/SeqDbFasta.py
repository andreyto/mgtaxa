"""A simple "database" of multi-FASTA files."""

from MGT.DirStore import *
from MGT.Taxa import *
from MGT.FastaIO import *
from MGT.FeatCommon import *
from MGT.FastaSplitters import *

class SeqDbFasta(DirKeyStore):
    """An interface to a collection of FASTA files.
    Each file is named with some ID, such as taxonomy ID.
    Methods for creating the DB grouped by taxonomy and streaming
    based on lists of IDs are provided.
    """

    objSfx = ".fasta.gz"
    objUncomprSfx = ".fasta"

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def setTaxaTree(self,taxaTree):
        self.taxaTree = taxaTree


    def loadIdList(self,objSfx=None):
        self.idList = self.listIds(objSfx=objSfx)

    def getIdList(self,objSfx=None):
        if not hasattr(self,"idList"):
            self.loadIdList(objSfx=objSfx)
        return self.idList

    def updateMetaDataById(self,id,force=False):
        meta = self.loadMetaDataById(id,{})
        if "seqLengths" not in meta or force:
            meta["seqLengths"] = self.computeSeqLengths(id)
        self.saveMetaDataById(id,meta)

    def finById(self,id):
        """Finalize creation of a new ID in SeqDb after using methods such as importByTaxa.
        This creates meta-data and possibly does other things."""
        #first compress, so that updateMetaDataById could find the file
        dataFileUncompr = self.getFilePathById(id,objSfx=self.objUncomprSfx)
        if os.path.exists(dataFileUncompr):
            compressFile(dataFileUncompr,self.getFilePathById(id))
        self.updateMetaDataById(id)

    def fastaReader(self,id):
        """Return FastaReader to read from the DB for a given ID.
        Caller is responsible for closing the reader."""
        return FastaReader(self.getFilePathById(id))

    def fastaReaders(self,ids):
        """Return iterator of FastaReader objects to read from the DB for a given set of IDs.
        Caller is responsible for closing the readers"""
        for id in ids:
            yield self.fastaReader(id)

    def fastaReaderChain(self,ids):
        """Return FastaReaderChain to read from the DB for a given set of IDs"""
        return FastaReaderChain(self.fastaReaders(ids))
    
    def computeSeqLengths(self,id):
        """Compute a numpy recarray with fields ("id","len") for a given DB ID
        directly from sequence data.
        Currently has to read all sequences for id - works through FastaReader"""
        return fastaLengths(self.getFilePathById(id))

    def seqLengths(self,id):
        """Return a numpy recarray with fields ("id","len") for a given DB ID.
        Loads from meta data that has to be pre-computed."""
        meta = self.loadMetaDataById(id)
        return meta["seqLengths"]
    
    def seqLengthSum(self,id):
        return self.seqLengths(id)["len"].sum()
    
    def seqCount(self,id):
        return len(self.seqLengths(id))

    def seqLengthsMany(self,ids=None):
        if ids is None:
            ids = self.iterIds()
        for id in ids:
            yield (id,self.seqLengths(id))

    def seqLengthsSumMany(self,ids=None):
        if ids is None:
            ids = self.iterIds()
        for id in ids:
            yield (id,self.seqLengthSum(id))
    
    def seqLengthsSumTotal(self,ids=None):
        return sum(self.seqLengthsSumMany(ids))
    
    def seqCountTotal(self,ids=None):
        if ids is None:
            ids = self.iterIds()
        return sum((self.seqCount(id) for id in ids))

    def cacheSeqCount(self):
        self.setAttr("seq_count",self.seqCountTotal())
    
    def fastaWriter(self,id,lineLen=None,mode="w",objSfx=None,compresslevel=6):
        """Return FastaWriter to write INTO the DB for a given ID.
        @param id ID of the record
        @param lineLen Length of FASTA lines
        @param mode String with the same semantics as 'mode' supplied to built-in open() method
        @param objSfx set to self.objUncomprSfx to write into uncompressed object during initial
        serial DB construction
        """
        return FastaWriter(out=self.getFilePathById(id,objSfx=objSfx),
                lineLen=lineLen,mode=mode,compresslevel=compresslevel)
    
    def fastaWriterUncompr(self,id,lineLen=None,mode="w",compresslevel=6):
        return self.fastaWriter(id,lineLen=lineLen,
                mode=mode,objSfx=self.objUncomprSfx,compresslevel=compresslevel)

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
                    #subsampled as a series of random fragments each.
                    #Look at preproc in FastaIO
                    #
                    #This gets a substring of a given length starting at a 
                    #random position within the full sequence
                    seq = SubSamplerRandomStart(int(len(seq)*lenRatio))(seq)
                if len(seq) > 0:
                    if fsout is None: 
                        fsout = FastaWriter(out,lineLen=reader.lineLen())
                    fsout.record(header=header,sequence=seq)
                    fsout.record(header=header,sequence=revCompl(seq))
                    #print "DEBUG:     id=%s seqLen=%s seqLenWritten=%s" % (id,seqLen,len(seq))
                    totalWritten += len(seq)
            reader.close()
        #print "DEBUG: id=%s maxLen=%s lenRatio=%s totalWritten=%s" % (id,maxLen,lenRatio,totalWritten)

