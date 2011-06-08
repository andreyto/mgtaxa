### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""I/O of FASTA formatted sequence files"""

__all__ = [ "FastaReader", "fastaReaderGzip", "FastaWriter", "splitFasta", "fastaLengths", "seqToLines" ]

from MGT.Common import *
from MGT.Util import openGzip, openCompressed
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
        if isinstance(infile,str):
            infile = openCompressed(infile,'r')
            self.ownInfile = True
        else:
            self.ownInfile = False
        self.infile = infile
        self.freshHdr = False
        self.maxLineLen = 0
        
    def records(self):
        infile = self.infile
        while True:
            if self.freshHdr:
                self.freshHdr = False
                yield self
                continue
            line = infile.next()
            # skip blank lines
            if line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                yield self
    
    def header(self):
        assert self.hdr.startswith('>')
        return self.hdr

    def getNCBI_Id(self):
        """Assume that header starts with '>gi|1234567|' and return the string id from second field."""
        return self.hdr.split('|',2)[1]
    
    def getNCBI_GI(self):
        return int(self.getNCBI_Id())

    def getSimpleId(self):
        """Assume that header starts with '>string_no_spaces ' and return that string."""
        return self.header().strip().split('>')[1].split(None,1)[0]

    def getId(self):
        """Attempt to extract and return an ID substring from the header.
        If the idExtractor method was set previously, use it.
        Otherwise, use getSimpleId()."""
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
            if line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                self.freshHdr = True
                return
            self.maxLineLen = max(self.maxLineLen,len(line)-1)
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

    def seqLen(self):
        n = 0
        for line in self.seqLines():
            n += len(line) - 1
            if not line.endswith("\n"):
                n += 1
        return n

    def lineLen(self):
        return self.maxLineLen

    def close(self):
        if self.ownInfile:
            self.infile.close()


def fastaReaderGzip(fileName):
    return FastaReader(openGzip(fileName,'r'))

def fastaLengths(inp):
    reader = FastaReader(inp)
    res = []
    for rec in reader.records():
        id = rec.getId()
        ln = rec.seqLen()
        res.append((id,ln))
    reader.close()
    return n.asarray(res,dtype=[("id",idDtype),("len","i8")])

def seqToLines(seq,lineLen=80):
    """Utility method that converts a string into line iterator.
    Endlines are inserted. The input string should not have endlines anywhere.
    This can be used in filters that emulate FastaReader record interface."""
    for x in range(0,len(seq),lineLen):
        yield seq[x:x+lineLen]+"\n"


class FastaWriter:
    
    def __init__(self,out,lineLen=None,mode="w"):
        if lineLen is None:
            lineLen = 1000
        if not hasattr(out,'write'):
            out = openCompressed(out,mode)
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


def splitFasta(inpFasta,outBase,maxChunkSize,sfxSep='_',lineLen=80):
    inpSeq = FastaReader(inpFasta)
    iChunk = 0
    chunkSize = 0
    def _open_new_out(iChunk):
        return openCompressed(outBase+'%s%04d'%(sfxSep,iChunk,),"w")
    out = _open_new_out(iChunk)
    print "Writing chunk %i of target size %i" % (iChunk,maxChunkSize)
    for rec in inpSeq.records():
        hdr = rec.header()
        out.write(hdr)
        lenSeq = 0
        #size rec.seqChunks() argument should be multiple of lineLen
        for chunk in rec.seqChunks(lineLen*1024):
            for x in range(0,len(chunk),lineLen):
                out.write(chunk[x:x+lineLen])
                out.write("\n")
            lenSeq += len(chunk)
            chunkSize += len(chunk)
        # we approximate the next seq length by the last one
        if chunkSize + lenSeq >= maxChunkSize:
            out.close()
            iChunk += 1
            chunkSize = 0
            out = _open_new_out(iChunk)
            print "Writing chunk %i of target size %i" % (iChunk,maxChunkSize)
    out.close()
    inpSeq.close()
    return iChunk

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

