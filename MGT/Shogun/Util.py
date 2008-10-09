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
    """Balance sample counts by to data['label'] and randomly shuffle entire set.
    @param data samples, must have data['label'] field
    @param maxCount if 0 - balance to the min class count, if <0 - only shuffle, else - to this count
    @param labTargets dict with optional per label maxCount values with the same semantics"""
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
            elif maxCount < 0:
                tcnt = cnt[lab]
            if tcnt > 0:
                dataLab = data[data['label'].astype('i4') == lab]
                dataLab = dataLab[random.sample(xrange(len(dataLab)),tcnt)]
                dataSel.append(dataLab)
    dataSel = numpy.concatenate(dataSel)
    #numpy permutation or shuffle do not work on arrays with 'O' datatypes
    dataSel = dataSel[nrnd.permutation(n.arange(len(dataSel),dtype=int))]
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

class LoadSeqPreprocShred:

    def __init__(self,sampLen,sampNum=0,sampOffset=0):
        self.sampLen = sampLen
        if isinstance(sampNum,int):
            sampNumConst = sampNum
            sampNum = lambda lab,seq,id: sampNumConst
        self.sampNum = sampNum
        self.sampOffset = sampOffset
        self.sampStride = sampLen + sampOffset
        assert sampLen > 0
        assert self.sampStride > 0

    def __call__(self,lab,seq,id):
        sampLen = self.sampLen
        sampStride = self.sampStride
        sampStartEnd = len(seq)-sampLen+1
        if sampStartEnd <= 0:
            return [],[],[]
        sampStarts = nrnd.permutation(n.arange(0,sampStartEnd,sampStride,dtype=int))
        sampNumVal = self.sampNum(lab,seq,id)
        if  sampNumVal > 0 and sampNumVal < len(sampStarts):
            sampStarts = sampStarts[:sampNumVal]
        sampSeq = [ seq[start:start+sampLen] for start in sampStarts ]
        return [lab]*len(sampSeq),sampSeq,[id]*len(sampSeq)

def loadSeqPreprocIdent(lab,seq,id):
    return ([lab], [seq], [id])

def loadSeqs(inpFile,preProc=loadSeqPreprocIdent,inpFileId=None):
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
        lab = float(lab)
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

def saveSeqs(data,outFile,outFileId=None):
    if outFileId is None:
        assert isinstance(outFile,str)
        outFileId = outFile+'.id'
    svmSaveId(data['id'],outFileId)
    if isinstance(outFile,str):
        outFile=openCompressed(outFile,'w')
        closeOut=True
    else:
        closeOut=False
    for rec in data:
        outFile.write("%s %s\n" % (rec['label'],rec['feature']))
    if closeOut:
        outFile.close()


