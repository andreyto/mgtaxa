### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""I/O of FASTA formatted sequence files"""

__all__ = [ "FastaReader", "fastaReaderGzip", "FastaWriter", "fastaLengths" ]

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
        if not hasattr(infile,"readline"):
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
            line = infile.readline()
            if not line:
                return
            # skip blank lines
            elif line.isspace():
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
        return self.hdr.split(None,1)[0][1:]

    def seqLines(self):
        infile = self.infile
        while True:
            line = infile.readline()
            if not line:
                break
            elif line.isspace():
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
        id = rec.header().split()[0][1:]
        ln = rec.seqLen()
        res.append((id,ln))
    reader.close()
    return n.asarray(res,dtype=[("id",idDtype),("len","i8")])

class FastaWriter:
    
    def __init__(self,out,lineLen=None):
        if lineLen is None:
            lineLen = 1000
        if not hasattr(out,'write'):
            out = openCompressed(out,'w')
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
        out = self.out
        if not header.startswith(">"):
            out.write(">")
        out.write(header)
        if not header.endswith("\n"):
            out.write("\n")
        if isinstance(sequence,str):
            s = sequence
        else:
            s = sequence.tostring()
        lineLen = self.lineLen
        for x in range(0,len(s),lineLen):
            out.write(s[x:x+lineLen])
            out.write("\n")

