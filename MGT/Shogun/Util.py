from MGT.Common import *
from MGT.Svm import *

class RevCompl:
    def __init__(self):
        trans = ['C']*256
        for (c,o) in zip('ATCG','TAGC'):
            trans[ord(c)] = o
        self.trans = ''.join(trans)

    def __call__(self,s):
        return s.translate(self.trans)[::-1]

def transDegen(seq):
    abet = 'ATCG'
    if not checkSaneAlphaHist(seq,nonDegenSymb=abet,minNonDegenRatio=0.9):
        print "Apparently, wrong alphabet for sequence: " + seq
    nAbet = len(abet)
    s = numpy.fromstring(seq.upper(),dtype='S1')
    for i in xrange(len(s)):
        if s[i] not in abet:
            #print i,s[i],"->",
            s[i] = abet[nrnd.randint(nAbet)]
            #print s[i]
    return s.tostring()

def randomSeq(lenSeq):
    return fromstring('ACTG',dtype='S1')[nrnd.randint(0,4,lenSeq)].tostring()

revCompl = RevCompl()

def alphaHist(seq):
    d = {}
    for c in seq:
        try:
            d[c] += 1
        except KeyError:
            d[c] = 1
    return d

def checkSaneAlphaHist(seq,nonDegenSymb,minNonDegenRatio=0.9):
    hist = numpy.bincount(numpy.fromstring(seq,dtype='b'))
    ind_nd = numpy.fromstring(nonDegenSymb,dtype='b')
    return float(hist[ind_nd[ind_nd<len(hist)]].sum())/len(seq) >= minNonDegenRatio

def checkSaneAlphaHistOld(seq,nonDegenSymb,minNonDegenRatio=0.9):
    h = alphaHist(seq)
    nNonD = sum([ h[c] for c in nonDegenSymb if c in h ])
    return float(nNonD)/len(seq) >= minNonDegenRatio

def balance(data,maxCount=0,labTargets={}):
    cnt = numpy.bincount(data['label'].astype('i4'))
    targCnt = cnt[cnt > 0].min()
    if maxCount > 0:
        targCnt = maxCount
    dataSel = []
    for lab in xrange(len(cnt)):
        if cnt[lab] > 0:
            tcnt = min(targCnt,cnt[lab])
            if lab in labTargets:
                labTarg = labTargets[lab]
                if labTarg >= 0:
                    tcnt = min(labTargets[lab],cnt[lab])
                else:
                    tcnt = cnt[lab]
            if tcnt > 0:
                dataLab = data[data['label'].astype('i4') == lab]
                dataLab = dataLab[random.sample(xrange(len(dataLab)),tcnt)]
                dataSel.append(dataLab)
    dataSel = numpy.concatenate(dataSel)
    return dataSel

def addRevComplRows(data):
    dataRC = data.copy()
    for rec in dataRC:
        rec['feature'] = revCompl(rec['feature'])
    return numpy.concatenate([data,dataRC])

def addRevComplCols(data):
    for rec in data:
        rec['feature'] = rec['feature']+revCompl(rec['feature'])
    return data

def applyRevCompl(data):
    for rec in data:
        rec['feature'] = revCompl(rec['feature'])
    return data

def splitStringFeat(data,sampLen):
    #dataRC = data.copy()
    applyToFeatData(data,lambda s: s[:sampLen])
    #applyToFeatData(dataRC,lambda s: s[sampLen:])
    return data
    #return numpy.concatenate([data,dataRC])

def applyToFeatData(data,func):
    for rec in data:
        rec['feature'] = func(rec['feature'])

def loadSeqs(inpFile,preProc=lambda lab,seq,id: ([lab], [seq], [id]),inpFileId=None):
    if inpFileId is None:
        assert isinstance(inpFile,str)
        inpFileId = inpFile+'.id'
    idsInp = svmLoadId(inpFileId)
    if isinstance(inpFile,str):
        inpFile=openCompressed(inpFile,'r')
        closeInp=True
    else:
        closeInp=False
    labs = []
    seqs = []
    ids  = []
    iLine = 0
    for line in inpFile:
        lab, seq = line.split(None,1)
        lab = int(lab)
        if seq[-1] == '\n':
            seq = seq[:-1]
        rec = preProc(lab,seq,idsInp[iLine])
        if rec is not None:
            labs.extend(rec[0])
            seqs.extend(rec[1])
            ids.extend(rec[2])
        iLine += 1

    label = numpy.asarray(labs,dtype='f8')
    feature = numpy.asarray(seqs,dtype='O')
    id = numpy.asarray(ids,dtype='O')
    data = numpy.rec.fromarrays((label,feature,id),names='label,feature,id')
    if closeInp:
        inpFile.close()
    return data

