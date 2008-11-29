#!/usr/bin/env python

from Bio import pairwise2
import numpy as n
import pdb

#A = "CTGCGGCCACGGTCGGTGCTGGTGGCCATAAAG"
#B = "CTGCGGCCACGGTCGGTGCTGGTGGCCATAAAG"
#B = "GACCCCGCCCTCTCCGACTTCGCTGGCGAGCCA"
#A = "ACCTGCGACCAGCACTGACGACCAGCCCAACCGGCGACTGTCACCGCAGGCCCGCGGCAGCAG"
A = "ACCTGCGACCAGCACTGACGACCAGCCCAACCGGCGACTGTCAGCAG"
#B = "TGGCGCGGGCGGCCATCGCTACCAGCGGCACGAGGGACCATCAACACCCGTTTGCGTCGGCATC"
B = "GCCTGCGAGCGGCGCCGACGATCAGCCCCGCCGGCAAGCGTCACCGCCCGTGCGTGGCCA"

print "Seq len: ", len(A), len(B)

def _pretty_print_align(align1, align2, score, begin, end):
    s = pairwise2.format_alignment(align1, align2, score, begin, end)
    a1 = n.fromstring(align1,dtype='S1')
    a2 = n.fromstring(align2,dtype='S1')
    print "Identity: %.2f Alignment length: %s" % (float(score)/len(align1)*100,len(align1))
    print "(a1 == a2).sum() = ", (a1 == a2).sum()
    print s,

def _pretty_print_all(aligns):
    aligns.sort()
    for align in aligns[:3]:
        _pretty_print_align(*align)

def _align_and_print(fn, *args, **keywds):
    print fn.__name__
    _pretty_print_all(fn(*args, **keywds))

a = pairwise2.align

print "### Test the function generation code"
fn = a.globalxx
print fn.__name__     # globalxx
print fn.__doc__      # globalxx(sequenceA, sequenceB) -> alignments
print a.localcd.__doc__  # localxx(...match_fn, openA, ...) -> alignments
try:
    a.blah
except:
    print "correctly failed"
    
print "#### Test penalize_end_gaps"
_align_and_print(a.globalxs, A, B, -0.8,-0.2,penalize_end_gaps=1) 

print "#### Let's start with a simple global alignment."
_align_and_print(a.globalxx, A, B)                 # 2 aligns, 3

print "#### Try a local alignment."
_align_and_print(a.localxs, A, B, -0.1, 0)         # 2 aligns, 1.9

print "#### Test match score, open penalty."
_align_and_print(a.globalms, A, B, 2.0, -1, -0.1, 0)    # 2 aligns, 1.9
_align_and_print(a.globalms, A, B, 1.5, 0, -0.1, 0)   # 2 aligns, 2.9
_align_and_print(a.globalxs, A, B, -0.1, 0)        # 1 align, 2.9
_align_and_print(a.globalms, A, B, 1, -2, -0.1, 0)  # 1 align, -0.1

print "#### Test the extend penalty."
_align_and_print(a.globalxs, A, B, -0.2, -0.5) # 1 align, 1.3
_align_and_print(a.globalxs, A, B, -0.2, -1.5) # 2 aligns, 0.6

print "#### Test penalize_extend_when_opening"
_align_and_print(a.globalxs, A, B, -0.2, -1.5,
                 penalize_extend_when_opening=1)       # 1 align, -1.2

print "#### Test penalize_end_gaps"
_align_and_print(a.globalxs, A, B, -0.2, -0.8,
                 penalize_end_gaps=0)                  # 3 aligns, 1

print "#### Test separate gap penalties"
_align_and_print(a.localxd, A, B, -0.3, 0, -0.8, 0)  # 2 aligns, 1.7
_align_and_print(a.localxd, A, B, -0.5, 0, -0.2, 0)  # 1 aligns, 1.8

#print "#### Test separate gap penalties, with extension.  Test align list"
#_align_and_print(a.localxd, list("GAAT"), list("GTCCT"), -0.1, 0, -0.1, -0.1,
#                 gap_char=["-"])                              # 3 aligns, 1.9

