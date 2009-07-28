"""CRISPR arrays"""

from MGT.Common import *

from Bio import pairwise2
from Bio.Blast import NCBIXML
import pdb

def crisprLengthRatioOk(arr,spacerToRepeatRatio=(0.6,2.5)):
    spLen = n.asarray([len(rec['spacer']) for rec in arr[:-1]])
    repLen = n.asarray([len(rec['repeat']) for rec in arr])
    avgRepLen = repLen.mean()
    if not n.logical_and(spLen >= avgRepLen*spacerToRepeatRatio[0],spLen <= avgRepLen*spacerToRepeatRatio[1]).all():
        return False
    return True

def crisprHasSimilarSpacers(arr,nLookAhead=3,maxIdentity=0.6,maxSimilarRatio=0.2):
    aligner = pairwise2.align.globalxs
    # Copy of comments from Bio.pairwise2 code:
    # The alignment functions take some undocumented keyword parameters:
    # - penalize_extend_when_opening: boolean
    #   Whether to count an extension penalty when opening a gap.  If
    #   false, a gap of 1 is only penalize an "open" penalty, otherwise it
    #   is penalized "open+extend".
    # - penalize_end_gaps: boolean
    #   Whether to count the gaps at the ends of an alignment.  By
    #   default, they are counted for global alignments but not for local
    #   ones.
    # - gap_char: string
    #   Which character to use as a gap character in the alignment
    #   returned.  By default, uses '-'.
    # - force_generic: boolean
    #   Always use the generic, non-cached, dynamic programming function.
    #   For debugging.
    # - score_only: boolean
    #   Only get the best score, don't recover any alignments.  The return
    #   value of the function is the score.
    # - one_alignment_only: boolean
    #   Only recover one alignment.
    spacers = arr[:-1]
    found = 0
    for i in xrange(len(spacers)):
        for j in xrange(i+1,min(i+1+nLookAhead,len(spacers))):
            alment = aligner(\
                    spacers[i]['spacer'],
                    spacers[j]['spacer'],
                    -0.8,-0.2,penalize_end_gaps=1,one_alignment_only=0)
            align1, align2, score, begin, end = alment[0]
            #if spacers[i]['pos'] == 8627:
            #    pdb.set_trace()
            a1 = n.fromstring(align1,dtype='S1')
            a2 = n.fromstring(align2,dtype='S1')
            nident = (a1 == a2).sum()
            ident = float(nident)/len(a1)
            if ident >= maxIdentity:
                found+=1
                break
    return found > maxSimilarRatio * len(spacers)


def parseSpacerBlast(inFile,minAlignLen=20,minBitScore=40,maxEvalue=1e-5,out=None):
    if out is None:
        out = sys.stdout
    inp = openCompressed(inFile,"r")
    blastRecs = NCBIXML.parse(inp)
    spcnt = defdict(int)
    for blastRec in blastRecs:
        for alignment in blastRec.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length >= minAlignLen \
                        and hsp.bits >= minBitScore \
                        and hsp.expect <= maxEvalue:
                #if hsp.align_length >= 20 and \
                #        abs(hsp.align_length - blastRec.query_letters) <= 2 and \
                #        hsp.align_length - hsp.identities <= 2 and \
                #        hsp.bits >= 40:
                    #print str(sorted(hsp.__dict__.items()))
                    #print str(sorted(blastRec.__dict__.items()))
                    #print str(sorted(alignment.__dict__.items()))
                    spcnt[blastRec.query] += 1
                    print >> out, 'Query: ', blastRec.query
                    print >> out, 'Hit: ', alignment.hit_def
                    print >> out, hsp.query
                    print >> out, hsp.match
                    print >> out, hsp.sbjct
                    print >> out,   'Query lenght: ', blastRec.query_letters,\
                            'Alignment length: ', hsp.align_length,\
                            'Bit-score: ', hsp.bits,\
                            'Query Coverage: ', round(float(hsp.align_length)/blastRec.query_letters*100),\
                            'Mismatches: ', (hsp.align_length - hsp.identities),\
                            'Identity: ', round(float(hsp.identities)/hsp.align_length*100),\
                            'E-value: %e' % hsp.expect
                    print >> out, "\n\n"
                    #pdb.set_trace()
    print >> out, "Number of unique spacers with conforming hits: %s" % len(spcnt)
    print >> out, sorted(spcnt.items())

