### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.FeatIOTxt import *
from MGT.FeatIOPkl import *

def SvmSparseFeatureWriter(out,outId=None,format=defFeatIOFormat):
    if format == "txt":
        return SvmSparseFeatureWriterTxt(out=out,outId=outId)
    elif format == "pkl":
        return  SvmFeatureWriterPkl(out=out,outId=outId)
    else:
        raise ValueError(format)

def SvmStringFeatureWriter(out,outId=None,format=defFeatIOFormat,lineLen=None):
    if format == "txt":
        return SvmStringFeatureWriterTxt(out=out,outId=outId)
    elif format == "pkl":
        return SvmFeatureWriterPkl(out=out,outId=outId)
    elif format == "fasta":
        return SvmFastaFeatureWriterTxt(out=out,lineLen=lineLen)
    else:
        raise ValueError(format)

def loadSeqsIdDef(inpFile):
    return svmLoadId(featIdFileNameDef(inpFile))

class LoadSeqPreprocIdFilter:
    """LoadSeq preprocessor that loads only records present in idToLab map and sets 'label' fields from idToLab"""

    def __init__(self,idToLab):
        self.idToLab = idToLab

    def __call__(self,lab,seq,id):
        if id in self.idToLab:
            return ([self.idToLab[id]], [seq], [id])
        else:
            return None
    
class IdMap:
    """Map of old IDs to new IDs created as a result of e.g. sequence shredding"""

    def __init__(self,records):
        """Constructor.
        @param records numpy recarray mapping "oldId" to "id" field, 1-to-N, so that "id" is unique"""
        self.records = records

    def getRecords(self):
        return self.records

    def getIdToRec(self):
        return dict( (rec["id"],rec) for rec in self.records )

    def selByOldId(self,ids):
        """Return a new IdMap object that contains only records within a given old ids sequence"""
        ids = set(ids)
        recs = self.getRecords()
        ind = n.asarray([ ind for (ind,id) in it.izip(it.count(),recs["oldId"]) if id in ids ],dtype=int)
        return self.__class__(records=recs[ind])

class LoadSeqPreprocShred:

    def __init__(self,sampLen,sampNum=-1,sampOffset=0,makeUniqueId=False,sortByStarts=False):
        """Constructor.
        @param sampLen length of each output fragment
        @param sampNum if < 0 - return all sampLen long fragments, if other number - that many,
        otherwise it should be f(lab,seq,id) and return the number of shreds for a given sequence.
        @param sampOffset optional offset between end of each shred and start of next (<0 for overlap)
        @param makeUniqueId if True, generate new UUIDs for shreds - getIdMap() can be used to get the mapping and coords
        @note Results will be returned as a reference to input sequence data, not a copy.
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
        self.sortByStarts = sortByStarts
        assert sampLen > 0
        assert self.sampStride > 0

    def __call__(self,lab,seq,id):
        sampLen = self.sampLen
        sampStride = self.sampStride
        sampStartEnd = len(seq)-sampLen+1
        if sampStartEnd <= 0:
            self.sampStarts = n.asarray([],dtype=int)
            return [],[],[]
        sampStarts = nrnd.permutation(n.arange(0,sampStartEnd,sampStride,dtype=int))
        sampNumVal = self.sampNum(lab,seq,id)
        if  sampNumVal > 0 and sampNumVal < len(sampStarts):
            sampStarts = sampStarts[:sampNumVal]
        elif sampNumVal == 0:
            self.sampStarts = n.asarray([],dtype=int)
            return [],[],[]
        if self.sortByStarts:
            sampStarts.sort()
        sampSeq = [ seq[start:start+sampLen] for start in sampStarts ]
        if self.makeUniqueId:
            ids = genId(n=len(sampSeq))
        else:
            ids = zeroId(n=len(sampSeq),val=id)
        oldIds = zeroId(n=len(sampSeq),val=id)
        self.sampCoord.append(n.rec.fromarrays([oldIds,ids,sampStarts,sampStarts+sampLen],names="oldId,id,sampStart,sampEnd"))
        self.sampStarts = sampStarts
        return [lab]*len(sampSeq),sampSeq,ids

    def getLastSampStarts(self):
        """Return sample starts generated during the last call to __call__()"""
        return self.sampStarts

    def getIdMap(self):
        return IdMap(records=n.concatenate(self.sampCoord))

def loadSeqPreprocIdent(lab,seq,id):
    return ([lab], [seq], [id])

def loadSeqPreprocDropFeat(lab,seq,id):
    return ([lab], [None], [id])

def loadSeqs(inpFile,preProc=loadSeqPreprocIdent,inpFileId=None,genMissingId=False,format=defFeatIOFormat):
    idFileName = None
    if inpFileId is None:
        assert isinstance(inpFile,str)
        inpFileId = featIdFileNameDef(inpFile)
        idFileName = inpFileId
        if not os.path.exists(inpFileId):
            inpFileId = None
    if inpFileId is not None:
        idsInp = svmLoadId(inpFileId)
        _genId = lambda i: idsInp[i]
    else:
        if genMissingId:
            _genId = lambda i: genId()
        else:
            raise ValueError("Feature ID file not found: %s" % (idFileName,))
    if format == "txt":
        inp = SvmStringFeatureReaderTxt(inpFile)
    elif format == "pkl":
        inp = SvmFeatureReaderPkl(inpFile)
    else:
        raise ValueError(format)
    labs = []
    seqs = []
    ids  = []
    iLine = 0
    for lab,seq in inp:
        rec = preProc(lab,seq,_genId(iLine))
        if rec is not None:
            labs.extend(rec[0])
            seqs.extend(rec[1])
            ids.extend(rec[2])
        iLine += 1
    data = MGTSparseData.makeEmpty(len(labs))
    data["label"] = labs
    data["feature"] = seqs
    data["id"] = ids
    # The following code did not work because 'feature' was promoted by 'asarray()' to 2D type 'O' array every time
    # when each record of 'seqs' has equal size.
    #label = n.asarray(labs,dtype=labelDtype)
    #feature = n.asarray(seqs,dtype='O')
    #id = n.asarray(ids,dtype=idDtype)
    #data = n.rec.fromarrays((label,feature,id),) #names='label,feature,id')
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
    label = n.asarray(labs,dtype=labelDtype)
    feature = n.asarray(seqs,dtype='O')
    id = n.asarray(ids,dtype=idDtype)
    return n.rec.fromarrays((label,feature,id),names='label,feature,id')

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

def loadSeqsMany(inpFiles,preProc=loadSeqPreprocIdent,inpFilesId=None,format=defFeatIOFormat):
    assert not isinstance(inpFiles,str) # must be sequence of strings
    if inpFilesId is None:
        inpFilesId = [None]*len(inpFiles)
    return n.concatenate([ loadSeqs(inpFile=inpFile,preProc=preProc,inpFileId=inpFileId,format=format) for \
            (inpFile,inpFileId) in zip(inpFiles,inpFilesId) ])


def makeDenseFeatDtype(label,feature,id):
    return n.dtype([("label",label.dtype),
                ("feature",feature.dtype,feature.shape[-1]),
                ("id",id.dtype)])

def makeDenseFeatDtypeFromShape(nAttr):
    return n.dtype([("label",label.dtype),
                ("feature",MGTDenseRealFeatures.defDtype,nAttr),
                ("id",id.dtype)])

def makeDenseFeature(label,feature,id):
    return n.rec.fromarrays((label,feature,id),
            dtype=makeDenseFeatDtype(label,feature,id))

def sparseToDenseSeqs(data):
    """Convert sparse representation of the 'feature' field into dense matrix one.
    @param data Numpy array with dtype [('label',any),('feature','O'),('id',idDtype)]
    where each 'feature' record is Numpy array (('ind',int),('val',float)).
    @return Numpy record array with dtype that has the same 'label' and 'id' fields,
    but 'feature' field is now a dense matrix."""
    feat = data['feature']
    #sparse ind starts from 1
    n_feat = max(rec['ind'].max() for rec in feat)
    m = n.zeros((len(feat),n_feat),dtype=feat.dtype["val"])
    for iRec in xrange(len(feat)):
        rec = feat[iRec]
        m[iRec,(rec['ind']-1,)] = rec['val']
    return makeDenseFeature(data['label'],m,data['id'])

def denseToSparseSeqs(data):
    """Convert dense representation of the 'feature' field into sparse SVM one.
    @param[in] data Numpy array with dtype [('label',any),('feature',1D array),('id',idDtype)]
    @return data Numpy array with dtype [('label',any),('feature','O'),('id','O')]
    where each 'feature' record is Numpy array (('ind',int),('val',float)).
    """
    newData = MGTSparseData.makeEmpty(len(data))
    for (newRec,rec) in it.izip(newData,data):
        newRec["id"] = rec["id"]
        newRec["label"] = rec["label"]
        newRec["feature"] = featVecDenseToSparse(rec["feature"])
    return newData

def loadSparseSeqs(inpFile,inpFileId=None,genMissingId=False,format=defFeatIOFormat):
    if format == "pkl":
        preProc = loadSeqPreprocIdent
    elif format == "txt":
        preProc = loadSeqPreprocParseSparse
    else:
        raise ValueError(format)
    return loadSeqs(inpFile=inpFile,preProc=preProc,inpFileId=inpFileId,genMissingId=genMissingId,format=format)

def loadSparseSeqsMany(inpFiles,idLab=None,format=defFeatIOFormat):
    if idLab is None:
        if format == "pkl":
            preProc = loadSeqPreprocIdent
        elif format == "txt":
            preProc = loadSeqPreprocParseSparse
        else:
            raise ValueError(format)
        return loadSeqsMany(inpFiles=inpFiles,preProc=preProc,format=format)
    else:
        data = loadSeqsMany(inpFiles=inpFiles,preProc=LoadSeqPreprocIdFilter(idToLab=idLab.getIdToLab()),format=format)
        if format == "pkl":
            return data
        elif format == "txt":
            return convSeqsToSparseInPlace(data)
        else:
            raise ValueError(format)

def convSeqsToSparseInPlace(data):
    return convFeatInPlace(data,preProc=loadSeqPreprocParseSparse)


def loadSparseSeqsAsDense(inpFile,inpFileId=None):
    return sparseToDenseSeqs(loadSparseSeqs(inpFile,inpFileId))

def saveSeqs(data,outFile,outFileId=None,format=defFeatIOFormat,lineLen=None):
    writer = SvmStringFeatureWriter(out=outFile,outId=outFileId,format=format,lineLen=lineLen)
    for rec in data:
        writer.write(label=rec["label"],feature=rec["feature"],id=rec["id"])
    if isinstance(outFile,str):
        writer.close()

def saveSparseSeqs(data,outFile,outFileId=None,format=defFeatIOFormat):
    writer = SvmSparseFeatureWriter(out=outFile,outId=outFileId,format=format)
    for rec in data:
        writer.write(label=rec["label"],feature=rec["feature"],id=rec["id"])
    if isinstance(outFile,str):
        writer.close()

def saveDenseSeqsAsSparse(data,outFile,outFileId=None,format=defFeatIOFormat):
    writer = SvmSparseFeatureWriter(out=outFile,outId=outFileId,format=format)
    for rec in data:
        writer.write(label=rec["label"],feature=featVecDenseToSparse(rec["feature"]),id=rec["id"])
    writer.close()

def saveDenseSeqsAsXyz(data,out,outFileId=None):
    writer = SvmDenseFeatureWriterXyz(out=out,nSamp=len(data),outId=outFileId)
    for rec in data:
        writer.write(label=rec["label"],feature=rec["feature"],id=rec["id"])
    writer.close()


