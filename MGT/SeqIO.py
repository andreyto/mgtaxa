from MGT.FastaIO import FasterReader

def pullNCBISeq(inSeqs,seqIds):
    giMap = dict( ( (s.gi,s) for s in seqIds ) )
    for inSeq in inSeqs:
        inpSeq = FastaReader(inSeq)
        for rec in inpSeq.records():
            gi = rec.getNCBI_GI()
            if gi in giMap:
                seqId = giMap[gi]
                id = seqId.id
                

