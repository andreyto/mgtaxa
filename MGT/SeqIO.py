### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.FastaIO import *
from MGT.Util import *
from MGT.Svm import *

def pullNCBISeqIds(inSeqs,seqIds,outSeq):
    """Depricated: use SeqImportApp"""
    writer = SvmStringFeatureWriterTxt(outSeq)
    for seqId in seqIds:
        seqId.seqFound = False
    giMap = dict( ( (s.gi,s) for s in seqIds ) )
    for inSeq in inSeqs:
        inpSeq = FastaReader(inSeq)
        for rec in inpSeq.records():
            gi = rec.getNCBI_GI()
            if gi in giMap:
                seqId = giMap[gi]
                id = seqId.id
                seq = rec.sequence()
                seqLenFull = len(seq)
                if seqId.maxLen >= 0:
                    seq = SubSamplerRandomStart(seqId.maxLen)(seq)
                writer.write(gi, seq, id)
                seqId.seqFound = True
                seqId.seqLenFull = seqLenFull
                seqId.seqLenSamp = len(seq)
        inpSeq.close()
    writer.close()


def taxidFromPhyFastaHeader(hdr):
    """Extract taxid from a FASTA header formatted for PHY database"""
    return int(hdr.split("taxid:",1)[1].split(" ",1)[0])

def joinFastaByTaxid(inpFastaFile,outFastaFile):
    reader = fastaReaderGzip(inpFastaFile)
    writer = openGzip(outFastaFile,'w')
    taxidLast = -1
    for rec in reader.records():
        taxid = taxidFromPhyFastaHeader(rec.header())
        if taxid != taxidLast:
            writer.write(">%s\n" % (taxid,))
            taxidLast = taxid
        for line in rec.seqLines():
            writer.write(line)
    reader.close()
    writer.close()


class FakeSequenceChecker:
    """Impelemnts assertions to perform on fake sequence generated when options.debugFakeSequence is set.
    Content of incoming fake sequence is like this:
    >gi|2695852|emb|Y13263.1|ABY13263 Acipenser baeri mRNA for immunoglobulin heavy chain, clone ScH 112
    x2695852x2695852x2695852x2695852x2695852x2695852x2695852x2695852x2695852x2695852
    x2695852x2695852x2695852x2695852x2695852x2695852x2695852
    It length can be slightly larger than the length of the original sequence, because we never wrap
    the GI numbers."""

    sepChar = 'x'

    def __init__(self):
        self.iTest = 0
        self.iMismatch = 0
        self.gi2taxa = None

    def loadGiTaxa(self):
        if self.gi2taxa is None:
            self.gi2taxa = loadObj(options.taxaPickled)

    def checkArrayGi(self,gi,seq):
        giStr = str(gi)
        sDtype = 'S%s' % len(giStr)
        words = seq.tostring().split(self.sepChar)
        if len(words) > 2:
            assert ( numpy.asarray(words[1:-1],dtype=sDtype) == giStr ).all()
        if len(words) > 1:
            assert giStr.startswith(words[-1])
        assert giStr.endswith(words[0])

    def checkArrayTaxid(self,taxid,seq,spacer=None):
        sepChar = self.sepChar
        gi2taxa = self.gi2taxa
        if spacer is None:
            seqStr = ( seq.tostring(), )
        else:
            seqStr = seq.tostring().split(spacer)
        for s in seqStr:
            words = s.split(sepChar)
            # We need at least one 'x12345x' to be sure
            # that we extracted the full GI.
            # Otherwise, we just let it go.
            if len(words) > 2:
                giStr = words[1]
                sDtype = 'S%s' % len(giStr)
                try:
                    assert ( numpy.asarray(words[1:-1],dtype=sDtype) == giStr ).all()
                    assert giStr.startswith(words[-1])
                    assert giStr.endswith(words[0])
                    assert gi2taxa[int(giStr)] == taxid
                except AssertionError:
                    giSeq = int(giStr)
                    taxidSeq = gi2taxa[giSeq]
                    print "Mismatch: giSeq = %s taxidSeq = %s taxid = %s words[:2] = %s" % \
                            (giSeq,taxidSeq,taxid,words[:2])
                    self.iMismatch += 1
                self.iTest += 1
                if self.iTest % 10000 == 0:
                    print "Mismatch %s out of %s" % (self.iMismatch,self.iTest)


g_fakeSequenceChecker = None

def getFakeSequenceChecker():
    """Return a singleton FakeSequenceChecker object.
    Use this function instead of constructing FakeSequenceChecker directly,
    so that we will not have several copies of a large gi2taxa index in memory."""
    global g_fakeSequenceChecker
    if g_fakeSequenceChecker is None:
        g_fakeSequenceChecker = FakeSequenceChecker()

    return g_fakeSequenceChecker
