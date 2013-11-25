"""Collection of general purpose utilitites for working with sequences"""

from MGT.Common import *
from MGT.Util import SymbolRunsCompressor

def transDegen(seq):
    """@todo this can be made much faster by computing a mask of degen seq positions and assigning random array"""
    abet = 'ATCG'
    if not checkSaneAlphaHist(seq,nonDegenSymb=abet,minNonDegenRatio=0.9):
        print "Apparently, wrong alphabet for sequence: " + seq.tostring()
    nAbet = len(abet)
    s = ArrStr.upper(seq)
    for i in xrange(len(s)):
        if s[i] not in abet:
            #print i,s[i],"->",
            s[i] = abet[nrnd.randint(nAbet)]
            #print s[i]
    return s

def randomSeq(lenSeq):
    return n.fromstring('ACTG',dtype='S1')[nrnd.randint(0,4,lenSeq)].view(n.chararray)


def alphaHist(seq):
    d = defdict(int)
    for c in seq:
        d[c] += 1
    return d

def checkSaneAlphaHist(seq,nonDegenSymb,minNonDegenRatio=0.9):
    if isinstance(seq,str):
        seq = n.fromstring(seq,dtype='S1')
    hist = n.bincount(n.asarray(seq.view('b')))
    ind_nd = n.fromstring(nonDegenSymb,dtype='b')
    return float(hist[ind_nd[ind_nd<len(hist)]].sum())/len(seq) >= minNonDegenRatio

def checkSaneAlphaHistOld(seq,nonDegenSymb,minNonDegenRatio=0.9):
    h = alphaHist(seq)
    nNonD = sum([ h[c] for c in nonDegenSymb if c in h ])
    return float(nNonD)/len(seq) >= minNonDegenRatio

