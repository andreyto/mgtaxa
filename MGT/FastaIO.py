### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""I/O of FASTA formatted sequence files"""

__all__ = [ "FastaReader", "FastaReaderChain", "fastaReaderGzip", "FastaWriter", 
        "splitFasta", "fastaLengths", "seqToLines", "shredFasta",
        "filterFastaByLength","fastaReaderFilterNucDegen",
        "estimateFastaLengthAndCountTotal" ]

from MGT.Common import *
from MGT.Util import openGzip, openCompressed, SymbolRunsCompressor, is_string
from MGT.SeqUtil import checkSaneAlphaHist
from MGT.UUID import *

from cStringIO import StringIO
import numpy

class FastaReader(object):
    """Class that supports an input iterator protocol for a FASTA file.
    Example that prints an exact copy of the input file:
    for rec in FastaReader(open('seq.fsa','r')).records():
        print rec.header(),
        for line in rec.seqLines():
            print line,
    Instead of rec.seqLines(), you can use the methods which post-process the
    raw sequences lines: seqChunks(), seqArrays(), sequence().
    """

    def __init__(self,infile):
        """Ctor.
        @param infile It can be either a string with a file name, or it can be an
        iterator that returns lines (each line should be terminated with a new line),
        e.g. a file object. An iterator can be a filter that reads another FastaReader
        object, performs transformations on the records and emits them as lines."""
        if is_string(infile):
            infile = openCompressed(infile,'r')
            self.ownInfile = True
        else:
            self.ownInfile = False
        self.infile = infile
        self.freshHdr = False
        self.maxLineLen = 0
        self.seqTotal = 0
        self.symbolsTotal = 0
        
    def records(self):
        infile = self.infile
        while True:
            if self.freshHdr:
                self.freshHdr = False
                yield self
                continue
            line = infile.next()
            self.symbolsTotal += len(line)
            # skip blank lines
            if line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                yield self
            else:
                #line can never be ""
                self.seqTotal += \
                        len(line) - 1 if line[-1].isspace() \
                        else len(line)

    
    def header(self):
        assert self.hdr.startswith('>')
        assert not self.freshHdr, "Detected invalid sequence of calls. Header "+\
                "already points to next record, but you did not advance your "+\
                "record iterator to it. You probably called "+\
                "rec.header(); rec.sequence(); rec.getId()."
        return self.hdr

    def getNCBI_Id(self):
        """Assume that header starts with '>gi|1234567|' and return the string id from second field."""
        return self.header().split('|',2)[1]
    
    def getNCBI_GI(self):
        return int(self.getNCBI_Id())

    def getSimpleId(self):
        """Assume that header starts with '>string_no_spaces ' and return that string."""
        return self.header().strip().split('>')[1].split(None,1)[0]

    def getId(self):
        """Attempt to extract and return an ID substring from the header.
        If the idExtractor method was set previously, use it.
        Otherwise, use getSimpleId().
        Note that this will work on whatever """
        if hasattr(self,"idExtractor"):
            return self.idExtractor(self.header())
        else:
            return self.getSimpleId()

    def setIdExtractor(self,idExtractor):
        """Assign a method to be used in future calls to getId().
        @param idExtractor A method that takes a header string and return sequences ID from it
        """
        self.idExtractor = idExtractor

    def seqLines(self):
        infile = self.infile
        while True:
            line = infile.next()
            self.symbolsTotal += len(line)
            if line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                self.freshHdr = True
                return
            lineLen = len(line) - 1 if line[-1].isspace() \
                    else len(line)
            self.maxLineLen = max(self.maxLineLen,lineLen)
            self.seqTotal += lineLen
            yield line

    def seqChunks(self,chunkSize):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
            if seq.tell() >= chunkSize:
                yield seq.getvalue()
                seq.close()
                seq = StringIO()
        if seq.tell() > 0:
            yield seq.getvalue()
        seq.close()

    def seqArrays(self,chunkSize):
        for s in self.seqChunks(chunkSize):
            yield numpy.fromstring(s,dtype='S1')

    def sequence(self,format='str'):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
        s = seq.getvalue()
        seq.close()
        if format == 'array':
            s = numpy.fromstring(s,dtype='S1')
        elif format != 'str':
            raise ValueError("Unknown format value: %s" % (format,))
        return s

    def seqLen(self,exclSymb=''):
        n = 0
        for line in self.seqLines():
            n += len(line.replace(exclSymb,'')) - 1
            if not line.endswith("\n"):
                n += 1
        return n

    def lineLen(self):
        return self.maxLineLen

    def seqLenLoaded(self):
        """Total number of sequence symbols loaded so far"""
        return self.seqTotal

    def symbolsLenLoaded(self):
        """Total number of symbols loaded so far.
        This includes headers, sequences and new lines"""
        return self.symbolsTotal

    def close(self):
        if self.ownInfile:
            self.infile.close()

class FastaReaderChain(object):
    """Class that chains multiple FastaReaders similar to itertools.chain().
    It emulates the interface of FastaReader to provide a drop-in replacement.
    """

    def __init__(self,fastaReaders):
        """Ctor.
        @param fastaReaders iterable of objects implementing FastaReader
        interface (could be FastaReaderChain objects too) or objects
        that can be passed to FastaReader constructor (file names or
        file objects).
        """
        self.fastaReaders = fastaReaders
        self.seqTotal = 0
        self.symbolsTotal = 0
        self.reader = None

    def records(self):
        #all other methods are provided by the returned record objects
        for reader in self.fastaReaders:
            if not hasattr(reader,"records"):
                reader = FastaReader(reader)
            self.reader = reader
            try:
                for rec in reader.records():
                    yield rec
            finally:
                self.seqTotal += reader.seqLenLoaded()
                self.symbolsTotal += reader.symbolsLenLoaded()
                reader.close()
                self.reader = None

    def close(self):
        #we already closed all readers through which we iterated
        #we should not try to close the remaining readers that the 
        #input iterator might yield - otherwise it would not be
        #possible for the client to stop iteration without
        #still generating all reader objects
        pass
    
    def seqLenLoaded(self):
        """Total number of sequence symbols loaded so far"""
        x = self.seqTotal
        if self.reader:
            x += self.reader.seqLenLoaded()
        return x

    def symbolsLenLoaded(self):
        """Total number of symbols loaded so far.
        This includes headers, sequences and new lines"""
        x = self.symbolsTotal
        if self.reader:
            x += self.reader.symbolsLenLoaded()
        return x


def fastaReaderGzip(fileName):
    return FastaReader(openGzip(fileName,'r'))

def fastaLengths(inp,exclSymb=''):
    reader = FastaReader(inp)
    res = []
    for rec in reader.records():
        id = rec.getId()
        ln = rec.seqLen(exclSymb=exclSymb)
        res.append((id,ln))
    reader.close()
    return n.asarray(res,dtype=[("id",idDtype),("len","i8")])

def estimateFastaLengthAndCountTotal(inpFileName):
    """Crude method to estimate the number of bases without reading the entire file.
    In a patholgical case (one file is one record) this will read the entire file."""
    symbolsToRead = 1024*1024
    with closing(FastaReader(inpFileName)) as inp:
        symb = 0
        nrec = 0
        for record in inp.records():
            nrec += 1
            symb = inp.symbolsLenLoaded()
            if symb >= symbolsToRead:
                break
        seq = inp.seqLenLoaded()
        fileSize = os.path.getsize(inpFileName)
        if getCompressionFormat(inpFileName):
            fileSize *= 3 #very crude
        if symb != 0:
            seqByTot = seq/float(symb)
        else:
            seqByTot = 1.
        if fileSize == 0:
            nSeqs = 0
        else:
            nSeqs = nrec * symb / fileSize
        #assumes that file in one byte encoded
        seqSize = fileSize * seqByTot
        return (seqSize,nSeqs)


def seqToLines(seq,lineLen=80):
    """Utility method that converts a string into line iterator.
    Endlines are inserted. The input string should not have endlines anywhere.
    This can be used in filters that emulate FastaReader record interface."""
    for x in range(0,len(seq),lineLen):
        yield seq[x:x+lineLen]+"\n"


class FastaWriter:
    
    def __init__(self,out,lineLen=None,mode="w",compresslevel=6):
        if lineLen is None:
            lineLen = 1000
        if not hasattr(out,'write'):
            out = openCompressed(out,mode,compresslevel=compresslevel)
            self.outClose = True
        else:
            self.outClose = False
        self.out = out
        self.lineLen = lineLen
    
    def __del__(self):
        if self.outClose:
            self.close()

    def close(self):
        self.out.close()

    def record(self,header,sequence):
        self.header(header)
        self.sequence(sequence)

    def header(self,header):
        out = self.out
        if not header.startswith(">"):
            out.write(">")
        out.write(header)
        if not header.endswith("\n"):
            out.write("\n")

    def sequence(self,sequence):
        out = self.out
        if isinstance(sequence,str):
            s = sequence
        else:
            s = sequence.tostring()
        lineLen = self.lineLen
        for x in range(0,len(s),lineLen):
            out.write(s[x:x+lineLen])
            out.write("\n")
    
    def seqLines(self,iterLines):
        out = self.out
        line = None
        for line in iterLines:
            out.write(line)
        if line and not line.endswith("\n"):
            out.write("\n")


def splitFasta(inpFasta,
        maxChunkSize,
        outBase=None,
        sfxSep=None,
        openWriter=None,
        lineLen=80,
        verbose=True):
    """Split input multi-FASTA file into multiple files of fixed size"""
    assert (outBase is None and sfxSep is None) or openWriter is None,\
            "Incompatible argument combination"
    if outBase is None:
        outBase = ""
    if sfxSep is None:
        sfxSep = '_'
    inpClose = False
    if hasattr(inpFasta,"records"):
        records = inpFasta.records()
    elif hasattr(inpFasta,"read") or is_string(inpFasta):
        inpFasta = FastaReader(inpFasta)
        inpClose = True
        records = inpFasta.records()
    else:
        records = inpFasta
    try:
        iChunk = 0
        chunkSize = 0
        out = None
        if openWriter is None:
            def _open_new_out(iChunk):
                return FastaWriter(outBase+'%s%04d'%(sfxSep,iChunk,),
                        lineLen=lineLen,mode="w")
        else:
            _open_new_out = openWriter
        out = _open_new_out(iChunk)
        try:
            if verbose:
                print "Writing chunk %i of target size %i" % (iChunk,maxChunkSize)
            for rec in records:
                hdr = rec.header()
                idSeq = rec.getId()
                seq = rec.sequence()
                out.record(hdr,seq)
                lenSeq = len(seq)
                yield (idSeq,lenSeq)
                chunkSize += lenSeq
                # we approximate the next seq length by the last one
                if chunkSize + lenSeq >= maxChunkSize:
                    out.close()
                    out = None
                    iChunk += 1
                    chunkSize = 0
                    out = _open_new_out(iChunk)
                    if verbose:
                        print "Writing chunk %i of target size %i" % (iChunk,maxChunkSize)
        finally:
            if out is not None:
                out.close()
    finally:
        if inpClose:
            inpFasta.close()

def shredFasta(inpFasta,outFasta,fragSize,fragCountRatio=1.,lineLen=80,outMode="w"):
    """Shred each record in multi-FASTA file into multiple records of fixed size"""
    from MGT.FeatIO import LoadSeqPreprocShred
    inpSeq = FastaReader(inpFasta)
    outSeq = FastaWriter(out=outFasta,lineLen=lineLen,mode=outMode)
    if fragCountRatio < 1.:
        sampNum = lambda lab,seq,id: int(len(seq)/(fragSize/2)*fragCountRatio) 
    else:
        sampNum = 0
    shredder = LoadSeqPreprocShred(sampLen=fragSize,sampNum=sampNum,
            sampOffset=-fragSize/2,makeUniqueId=False,
            sortByStarts=True)
    for rec in inpSeq.records():
        hdr = rec.header()
        id = rec.getId()
        seq = rec.sequence()
        labFr,seqFr,idFr = shredder(0,seq,id)
        startsFr = shredder.getLastSampStarts()
        for iF,sF in enumerate(seqFr):
            stF = startsFr[iF]
            idF = "%s_%s_%s-%s" % (id,iF,stF,stF+fragSize)
            outSeq.record(header=idF,sequence=sF)
    inpSeq.close()
    outSeq.close()

def fastaWriteOnce(records,out,lineLen=None,mode="w"):
    """Open out file, write several records at once as FASTA, and close the file.
    @param records Iterable of tuples (header,sequence)
    @param out Output file
    @param lineLen Line length for FASTA output
    @return Number of records written
    """
    w = FastaWriter(out=out,lineLen=lineLen,mode=mode)
    for (iRec,(hdr,seq)) in enumerate(records):
        w.record(header=hdr,sequence=seq)
    w.close()
    return (iRec+1)

def filterFastaByLength(inp,out,minLen=None,maxLen=None,lineLen=None,mode="w"):
    if minLen is None:
        minLen = 0
    if maxLen is None:
        maxLen = sys.maxint
    rdClose = False
    rd = inp
    wrClose = False
    wr = out
    try:
        if not hasattr(inp,"records"):
            rd = FastaReader(inp)
            rdClose = True
        try:
            if not hasattr(out,"record"):
                wr = FastaWriter(out=out,lineLen=lineLen,mode=mode)
                wrClose = True
            for rec in rd.records():
                hdr = rec.header()
                seq = rec.sequence()
                lenSeq = len(seq)
                if lenSeq >= minLen and lenSeq < maxLen:
                    wr.record(hdr,seq)
                    yield dict(rec=rec,length=lenSeq)
        finally:
            if wrClose:
                wr.close()
    finally:
        if rdClose:
            rd.close()

def fastaReaderFilterNucDegen(fastaReader,extraFilter=None,minNonDegenRatio=0.98,doWarn=True):
    if extraFilter is None:
        extraFilter = lambda hdr,seq: (hdr,seq)
    compr = SymbolRunsCompressor(sym="Nn",minLen=1)
    nonDegenSymb = "ATCGatcg"
    def line_gen():
        for rec in fastaReader.records():
            hdr = rec.header()
            seq = compr(rec.sequence())
            if not checkSaneAlphaHist(seq,nonDegenSymb,
                    minNonDegenRatio=minNonDegenRatio):
                if doWarn:
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
