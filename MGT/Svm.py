from MGT.Common import *
from MGT.Taxa import *
from MGT import Kmers
from itertools import izip

class SvmSparseFeatureWriterTxt(Kmers.SvmSparseFeatureWriterTxt):

    def __init__(self,out,outId=None):
        Kmers.SvmSparseFeatureWriterTxt.__init__(self,out)
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(out+'.id','w', buffering=1024*1024)
        self.outId = outId

    def write(self,label,feature,id=None):
        Kmers.SvmSparseFeatureWriterTxt.write(self,label,feature['values'],feature['indices'])
        if id is None:
            id = label
        self.outId.write("%s\n" % id)

    def close(self):
        Kmers.SvmSparseFeatureWriterTxt.close(self)
        self.outId.close()

class SvmFastaFeatureWriterTxt:
    
    def __init__(self,out):
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        
    def close(self):
        self.out.close()

    def write(self,label,feature):
        self.out.write(">%d\n" % label)
        self.out.write(feature.tostring())
        self.out.write("\n")
        self.nOut += 1

    def numRec(self):
        return self.nOut

def svmSaveId(ids,out):
    if not hasattr(out,'write'):
        out = openCompressed(out,'w')
        ownOut = True
    else:
        ownOut = False
    s = ''.join([ "%s\n" % x for x in ids ])
    out.write(s)
    if ownOut:
        out.close()

def svmLoadId(inp,dtype=idDtype):
    if not hasattr(inp,'read'):
        inp = openCompressed(inp,'r')
        ownInp = True
    else:
        ownInp = True
    x = n.asarray([ s.strip() for s in inp ],dtype=dtype)
    if ownInp:
        inp.close()
    return x
    

class SvmStringFeatureWriterTxt:
    
    def __init__(self,out,outId=None):
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(out+'.id','w', buffering=1024*1024)
        self.outId = outId
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        
    def close(self):
        self.out.close()
        self.outId.close()

    def write(self,label,feature,id=None):
        self.out.write("%d " % label)
        if isinstance(feature,str):
            self.out.write(feature)
        else:
            self.out.write(feature.tostring())
        self.out.write("\n")
        if id is None:
            id = label
        self.outId.write("%s\n" % id)
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmDenseFeatureWriterTxt:
    
    def __init__(self,out,outId=None):
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(out+'.id','w', buffering=1024*1024)
        self.outId = outId
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = None
        
    def close(self):
        self.out.close()
        self.outId.close()

    def write(self,label,feature,id=None):
        if self.formatStr is None:
            self.formatStr = ' '.join( ("%d:%%g" % ind for ind in xrange(1,len(feature)+1)) ) + '\n'
        self.out.write("%d " % (label,))
        self.out.write(self.formatStr % tuple(feature))
        if id is None:
            id = label
        self.outId.write("%s\n" % id)
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmDenseFeatureWriterCsv:
    
    def __init__(self,out,writeHeader=True,nFeat=None):
        assert not (writeHeader and nFeat is None)
        self.nFeat = nFeat
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = None
        if writeHeader:
            out.write("id,label,"+",".join([ "f%i" % iFeat for iFeat in xrange(nFeat) ])+"\n")
        
    def close(self):
        self.out.close()

    def write(self,label,feature,id=None):
        if self.formatStr is None:
            self.formatStr = ",%g"*len(feature)+'\n'
        if id is None:
            id = label
        self.out.write("%s,%d" % (id,label))
        self.out.write(self.formatStr % tuple(feature))
        self.nOut += 1

    def numRec(self):
        return self.nOut

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

def isUniqueId(data):
    return len(n.unique1d(data["id"])) == len(data["id"])

class IdMap:
    """Map of old IDs to new IDs created as a result of e.g. convFeat call"""

    def __init__(self,data):
        self.data = data

class IdLabels:
    """Maps unique sample IDs into classification labels (many-to-one relationship) and their splits"""

    def __init__(self,records=None,fileName=None,initMaps=True):
        """Ctor.
        @param records Numpy record array with fields id,label,split.
        @param fileName file created by save() call
        @pre Either records or fileName must be None. id field is unique, e.g. obtained with UUID module"""
        assert records is None or fileName is None
        self.data = Struct()
        if records is not None:
            self.data.records = n.unique1d(records)
            if initMaps:
                self.initMaps()
                self.check()
        elif fileName is not None:
            self.load(fileName,initMaps=initMaps)

    def initMaps(self):
        records = self.getRecords()
        self.idToRec = dict( (rec["id"],rec) for rec in records )
        self.labToRec = groupRecArray(records,keyField="label")
        self.splitToRec = groupRecArray(records,keyField="split")
        self.initLabToNameMap()

    def initLabToNameMap(self):
        if not hasattr(self.data,"labNames"):
            self.labToName = dict( ( (lab,lab) for lab in self.labToRec ) )

    def check(self):
        records = self.getRecords()
        idRecs = groupRecArray(records,keyField="id")
        nonUniqueIds = dict( ( (key,val) for (key,val) in idRecs.iteritems() if len(val) > 1 ) )
        assert len(nonUniqueIds) == 0
        
    def getRecords(self):
        return self.data.records

    def getIdToRec(self):
        return self.idToRec

    def getIdToLab(self):
        return dict( ( (rec["id"],rec["label"]) for rec in self.getRecords() ) )

    def getLabToRec(self):
        return self.labToRec

    def getSplitToRec(self):
        return self.splitToRec

    def getSplits(self):
        splits = {}
        splitToRec = self.getSplitToRec()
        for key,rec in splitToRec.iteritems():
            splits[key] = IdLabels(records=rec)
        return splits

    def getLabToName(self):
        return self.labToName

    def getLabNames(self):
        return self.data.labNames

    def setLabNames(self,labNames):
        self.data.labNames = labNames
        self.initLabToNameMap()

    def save(self,fileName):
        dumpObj(self.data,fileName)

    def load(self,fileName,initMaps=True):
        self.data = loadObj(fileName)
        if initMaps:
            self.initMaps()

    def union(self,other):
        """Return a union of self and other(s).
        @param other either IdLabels object or a sequence of such objects.
        @pre if there are identical ID values, the corresponding records are fully identical too."""
        srec = self.getRecords()
        if isinstance(other,IdLabels):
            orec = other.getRecords()
            nrec = n.concatenate([srec,orec])
        else:
            #other must be a sequence of IdLabels objects
            orec = [ x.getRecords() for x in other ]
            orec.append(srec)
            nrec = n.concatenate(orec)
        return self.__class__(records=nrec)

    def selDataInd(self,data):
        """Return array of row indices for data rows that are present in self"""
        idToRec = self.getIdToRec()
        return n.asarray([ ind for (ind,id) in it.izip(it.count(),data["id"]) if id in idToRec ],dtype='i4')

    def selData(self,data,setLab=True):
        d = data[self.selDataInd(data)]
        if setLab:
            self.setDataLab(d)
        return d

    def setDataLab(self,data):
        idToRec = self.getIdToRec()
        for rec in data:
            rec["label"] = idToRec[rec["id"]]["label"]


def unionIdLabels(idLabs):
    if len(idLabs) == 1:
        return idLabs[0]
    elif len(idLabs) > 1:
        first,rest = idLabs[0],idLabs[1:]
        return first.union(rest)
    else:
        raise ValueError("Empty sequence is not allowed as argument to unionIdLabels")


def loadIdLabelsMany(fileNames):
    return unionIdLabels([ IdLabels(fileName=fileName) for fileName in fileNames ])
        

def saveIdLabelRecords(records,fileName):
    """Save records that can be used to build an IdLabels object as such an object.
    This avoids spending time on building IdLabels internal maps."""
    idLabs = IdLabels(records=records,initMaps=False)
    idLabs.save(fileName)


class LoadSeqPreprocIdFilter:

    def __init__(self,idToLab):
        self.idToLab = idToLab

    def __call__(self,lab,seq,id):
        if id in self.idToLab:
            return ([self.idToLab[id]], [seq], [id])
        else:
            return None
    

class LoadSeqPreprocShred:

    def __init__(self,sampLen,sampNum=0,sampOffset=0,makeUniqueId=False):
        """Constructor.
        @param sampLen length of each output fragment
        @param sampNum if <=0 - return all sampLen long fragments, if other number - that many,
        otherwise it should be f(lab,seq,id) and return the number of shreds for a given sequence.
        @param sampOffset optional offset between end of each shred and start of next
        @param makeUniqueId if True, generate new UUIDs for shreds - getIdMap() can be used to get the mapping and coords
        """
        self.sampLen = sampLen
        if isinstance(sampNum,int):
            sampNumConst = sampNum
            sampNum = lambda lab,seq,id: sampNumConst
        self.sampNum = sampNum
        self.sampOffset = sampOffset
        self.sampStride = sampLen + sampOffset
        self.sampCoord = []
        self.makeUniqueId = makeUniqueId
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
        if self.makeUniqueId:
            ids = genId(n=len(sampSeq))
        else:
            ids = zeroId(n=len(sampSeq),val=id)
        oldIds = zeroId(n=len(sampSeq),val=id)
        self.sampCoord.append(n.rec.fromarrays([oldIds,ids,sampStarts],names="oldId,id,sampStart"))
        return [lab]*len(sampSeq),sampSeq,ids

    def getIdMap(self):
        return IdMap(data=n.concatenate(self.sampCoord))

def loadSeqPreprocIdent(lab,seq,id):
    return ([lab], [seq], [id])

def loadSeqPreprocDropFeat(lab,seq,id):
    return ([lab], [None], [id])

def loadSeqPreprocParseSparseOrig(lab,seq,id):
    feat = n.fromiter(( (int(ind),float(val)) 
        for ind,val in (entry.split(':') 
            for entry in seq.split()) ),dtype=[("ind","i4"),("val","f8")])
    return ([lab], [feat], [id])

_parseSparseTransTable = string.maketrans(':',' ')
#This one is about 70% faster than versions 1 & 2 below and 5x faster than the Orig
def loadSeqPreprocParseSparse(lab,seq,id):
    feat = n.fromstring(seq.translate(_parseSparseTransTable),sep=' ',dtype='f8')
    feat = n.rec.fromarrays([feat[::2],feat[1::2]],dtype=[("ind","i4"),("val","f8")])
    return ([lab], [feat], [id])

def loadSeqPreprocParseSparse1(lab,seq,id):
    feat = n.fromstring(seq.translate(_parseSparseTransTable),sep=' ',dtype='f8')\
            .view([("ind","f8"),("val","f8")]).astype([("ind","i4"),("val","f8")])
    return ([lab], [feat], [id])

def loadSeqPreprocParseSparse2(lab,seq,id):
    feat = n.rec.fromrecords([ entry.split(':') for entry in seq.split() ],dtype=[("ind","i4"),("val","f8")])
    return ([lab], [feat], [id])

def featIdFileNameDef(featFile):
    return featFile+'.id'

def loadSeqs(inpFile,preProc=loadSeqPreprocIdent,inpFileId=None,genMissingId=False):
    if inpFileId is None:
        assert isinstance(inpFile,str)
        inpFileId = featIdFileNameDef(inpFile)
        if not os.path.exists(inpFileId):
            inpFileId = None
    if inpFileId is not None:
        idsInp = svmLoadId(inpFileId)
        _genId = lambda i: idsInp[i]
    else:
        if genMissingId:
            _genId = lambda i: genId()
        else:
            raise ValueError("Feature ID file not found")
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
        rec = preProc(lab,seq,_genId(iLine))
        if rec is not None:
            labs.extend(rec[0])
            seqs.extend(rec[1])
            ids.extend(rec[2])
        iLine += 1
    label = numpy.asarray(labs,dtype='f8')
    feature = numpy.asarray(seqs,dtype='O')
    id = numpy.asarray(ids,dtype=idDtype)
    data = numpy.rec.fromarrays((label,feature,id),names='label,feature,id')
    if closeInp:
        inpFile.close()
    return data

def convFeat(data,preProc):
    labs = []
    seqs = []
    ids  = []
    for oldRec in data:
        rec = preProc(lab=oldRec["label"],seq=oldRec["feature"],id=oldRec["id"])
        if rec is not None:
            labs.extend(rec[0])
            seqs.extend(rec[1])
            ids.extend(rec[2])
    label = numpy.asarray(labs,dtype='f8')
    feature = numpy.asarray(seqs,dtype='O')
    id = numpy.asarray(ids,dtype=idDtype)
    return numpy.rec.fromarrays((label,feature,id),names='label,feature,id')

def convFeatInPlace(data,preProc):
    """Convert feature array in place.
    @param data feature array, will have its content destroyed
    @param preProc pre-processor object, that returns zero or one feature records on each call
    @ret resulting data array"""
    iPutRec = 0
    for oldRec in data:
        rec = preProc(lab=oldRec["label"],seq=oldRec["feature"],id=oldRec["id"])
        if rec is not None:
            assert len(rec[0]) == 1 and len(rec[1]) == 1 and len(rec[2]) == 1
            putRec = data[iPutRec]
            putRec["label"] = rec[0][0]
            putRec["feature"] = rec[1][0]
            putRec["id"] = rec[2][0]
            iPutRec += 1
    data.resize(iPutRec)
    return data

def loadSeqsMany(inpFiles,preProc=loadSeqPreprocIdent,inpFilesId=None):
    if inpFilesId is None:
        inpFilesId = [None]*len(inpFiles)
    return n.concatenate([ loadSeqs(inpFile=inpFile,preProc=preProc,inpFileId=inpFileId) for \
            (inpFile,inpFileId) in zip(inpFiles,inpFilesId) ])

def sparseToDenseSeqs(data):
    """Convert sparse representation of the 'feature' field into dense matrix one.
    @param data Numpy array with dtype [('label',any),('feature','O'),('id','O')]
    where each 'feature' record is Numpy array (('ind',int),('val',float)).
    @return Numpy record array with dtype that has the same 'label' and 'id' fields,
    but 'feature' field is now a dense matrix."""
    feat = data['feature']
    #sparse ind starts from 1
    n_feat = max(rec['ind'].max() for rec in feat)
    m = n.zeros((len(feat),n_feat),dtype='f8')
    for iRec in xrange(len(feat)):
        rec = feat[iRec]
        m[iRec,(rec['ind']-1,)] = rec['val']
    return n.rec.fromarrays((data['label'],m,data['id']),
            dtype=[("label",data.dtype["label"]),
                ("feature",m.dtype,m.shape[-1]),
                ("id",data.dtype["id"])])

def loadSparseSeqs(inpFile,inpFileId=None):
    return loadSeqs(inpFile=inpFile,preProc=loadSeqPreprocParseSparse,inpFileId=inpFileId)

def loadSparseSeqsMany(inpFiles,idLab=None):
    if idLab is None:
        return loadSeqsMany(inpFiles=inpFiles,preProc=loadSeqPreprocParseSparse)
    else:
        data = loadSeqsMany(inpFiles=inpFiles,preProc=LoadSeqPreprocIdFilter(idToLab=idLab.getIdToLab()))
        return convSeqsToSparseInPlace(data)

def convSeqsToSparseInPlace(data):
    return convFeatInPlace(data,preProc=loadSeqPreprocParseSparse)


def loadSparseSeqsAsDense(inpFile,inpFileId=None):
    return sparseToDenseSeqs(loadSparseSeqs(inpFile,inpFileId))

def saveSeqs(data,outFile,outFileId=None):
    if outFileId is None:
        assert isinstance(outFile,str)
        outFileId = featIdFileNameDef(outFile)
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


class LabelRenum:
    """Renumber arbitrary type labels to consequitive numbers"""

    def __init__(self,labels):
        labels = set(labels)
        labrenum = [ (lab,renum) for (lab,renum) in it.izip(sorted(labels),it.count()) ]
        self.oldToNew = dict(labrenum)
        self.newToOld = dict( ( (renum,lab) for (lab,renum) in labrenum ) )

    def toNew(self):
        return self.oldToNew

    def toOld(self):
        return self.newToOld
        
    def setDataLabels(self,data,labelsOld):
        oldToNew = self.toNew()
        for (rec,labOld) in it.izip(data,labelsOld):
            rec['label'] = oldToNew[labOld]

def setLabelsFromIds(data):
    labRenum = LabelRenum(data["id"])
    labRenum.setDataLabels(data,data["id"])
    return labRenum

class SVMLibLinear:
    
    def __init__(self,workDir):
        self.workDir = workDir
        self.modelFileRel = "model.svm"
        makedir(workDir)
        self.binDir = os.environ['SVM_LIBLINEAR_BIN']
        self.binTrain = os.path.join(self.binDir,'train')
        self.binPredict = os.path.join(self.binDir,'predict')
        
    def trainCmd(self,trainFile):
        #cmd = [self.binTrain,"-c","4","-e","0.1","-s","3",trainFile,self.modelFileRel]
        cmd = ["(zcat %s; echo '#'; zcat %s) |" % (trainFile,trainFile),
                self.binTrain,"-c","0.01","-e","1","-s","3","-",self.modelFileRel]
        return cmd
        #run(cmd, cwd=self.workDir)

        

class SVMLib:
    def __init__(self):
        import svm
        import cross_validation
        self.svm = svm
        self.cross_validation = cross_validation
        self.param = self.svm.svm_parameter(svm_type = self.svm.C_SVC, 
                                            kernel_type = self.svm.LINEAR, #self.svm.RBF,
                                            C = 2, gamma = 2, 
                                            cache_size=2000)
    
    def train(self,inpPath,testMode=False,testRecNum=10000):
        modelFile = inpPath + '.svm'
        kmersInp = KmerBinReader(inpPath)
        data = kmersInp.readAll()
        if testMode:
            data = numpy.concatenate((data[:testRecNum/2],data[-testRecNum/2:-1]))
        kmersInp.close()
        labels = data['taxid'].astype(numpy.float64)
        vals = data['vals'].astype(numpy.float64)
        print "labels: ", labels.shape, labels.dtype
        print "vectors: ", vals.shape, vals.dtype
        action = 'crossVal'
        if action == 'crossVal':
            print "Starting cross-validation with %s samples of rank %s" % (len(labels),len(vals[0]))
            self.cross_validation.do_cross_validation(vals,labels,self.param,5)
            print "Finished cross-validation"
        elif action == 'train':
            prob = self.svm.svm_problem(labels,vals)
            print "Starting training with %s samples of rank %s" % (len(labels),len(vals[0]))
            model = self.svm.svm_model(prob, self.param)
            print "Finished training"
            model.save(modelFile)
        else:
            raise ValueError("Unknown 'action' value: "+action)

class SVMMulticlass:
    
    def __init__(self,workDir='/export/tmp'):
        self.workDir = workDir
        self.sampleFileRel = "samples.svm"
        self.modelFileRel = "model.svm"
        makedir(workDir)
        
    def train(self,inpFile,testMode=False,testRecNum=10000):
        samples = open(os.path.join(self.workDir,self.sampleFileRel),'w', buffering=1024*1024)
        iRec = 0
        freqs = KmerFreqs(inpFile)    
        for rec in freqs.readValues():
            rec.write(samples)
            #samples.write(rec.asStr())
            samples.write("\n")
            iRec = iRec + 1
            if testMode and iRec >= testRecNum:
                break
        freqs.close()
        #important to flush, otherwise training program does not see the last records:
        samples.close()
        print "Starting training with %s samples of rank %s" % (iRec,len(rec.vals))        
        run(["svm_multiclass_learn","-c","1","-m","2000",self.sampleFileRel,self.modelFileRel], cwd=self.workDir)
        print "Finished training"



