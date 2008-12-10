from MGT.FastaIO import FastaReader
from MGT.Util import *
from MGT.Svm import *

def pullNCBISeq(inSeqs,seqIds,outSeq):
    writer = SvmStringFeatureWriterTxt(outSeq)
    giMap = dict( ( (s.gi,s) for s in seqIds ) )
    for inSeq in inSeqs:
        inpSeq = FastaReader(inSeq)
        for rec in inpSeq.records():
            gi = rec.getNCBI_GI()
            if gi in giMap:
                seqId = giMap[gi]
                id = seqId.id
                seq = rec.sequence()
                if seqId.maxLen >= 0:
                    seq = SubSamplerRandomStart(seqId.maxLen)(seq)
                writer.write(gi, seq, id)
        inpSeq.close()
    writer.close()

