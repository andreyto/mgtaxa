### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Taxa import *
from MGT.FeatIO import *
from MGT.FeatCommon import *


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
    hist = n.bincount(n.asarray(seq.view('b')))
    ind_nd = n.fromstring(nonDegenSymb,dtype='b')
    return float(hist[ind_nd[ind_nd<len(hist)]].sum())/len(seq) >= minNonDegenRatio

def checkSaneAlphaHistOld(seq,nonDegenSymb,minNonDegenRatio=0.9):
    h = alphaHist(seq)
    nNonD = sum([ h[c] for c in nonDegenSymb if c in h ])
    return float(nNonD)/len(seq) >= minNonDegenRatio

def balance(data,maxCount=0,labTargets={},rndSeed=None):
    """Balance sample counts by to data['label'] and randomly shuffle entire set.
    @param data samples, must have data['label'] field
    @param maxCount if 0 - balance to the min class count, if <0 - only shuffle, else - to this count, 
    if "median" - to median count
    @param labTargets dict with optional per label maxCount values with the same semantics"""
    # we want to use the global random module-level generator if a seed has not been provided to
    # this function, so that the globally set seed (if any) can still be in effect
    if rndSeed is not None:
        rngen = random.Random() #class instance
        rngen.seed(rndSeed)
    else:
        rngen = random #module, interface is same as class Random
    cnt = n.bincount(data['label'].astype('i4'))
    if isinstance(maxCount,str):
        if maxCount == "median":
            maxCount = int(round(n.median(cnt[cnt>0])))
        else:
            raise ValueError(maxCount)
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
                dataLab = dataLab[rngen.sample(xrange(len(dataLab)),tcnt)]
                dataSel.append(dataLab)
    dataSel = n.concatenate(dataSel)
    #we use already seeded generator from 'random' instead of possibly more efficient numpy.random
    ind = n.arange(len(dataSel),dtype=int)
    rngen.shuffle(ind)
    ###if we just used random.shuffle(dataSel), we would get dataSel sliently screwed up (some rows doubled, some gone) -
    ### seems to be a problem with swapping the numpy recarray rows in random.shuffl implementation.
    #rngen.shuffle(dataSel)
    dataSel = dataSel[ind]
    #this will have to be commented out if we ever use on the data array rather than IdLabels object.
    assert len(dataSel) == len(set(dataSel["id"]))
    ###numpy permutation or shuffle do not work on arrays with 'O' datatypes
    #dataSel = dataSel[nrnd.permutation(n.arange(len(dataSel),dtype=int))]
    return dataSel

def addRevComplRows(data):
    dataRC = data.copy()
    for rec in dataRC:
        rec['feature'] = revCompl(rec['feature'])
    return n.concatenate((data,dataRC))

def addRevComplCols(data):
    for rec in data:
        rec['feature'] = n.concatenate((rec['feature'],revCompl(rec['feature'])))
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
    #return n.concatenate([data,dataRC])

def applyToFeatData(data,func):
    for rec in data:
        rec['feature'] = func(rec['feature'])

def isUniqueId(data):
    return len(n.unique1d(data["id"])) == len(data["id"])


class IdLabels:
    """Maps unique sample IDs into classification labels (many-to-one relationship) and their splits"""

    _defRecDtype = n.dtype([("id",idDtype),("label",labelDtype),("split",splitDtype)])

    def __init__(self,records=None,fileName=None,initMaps=True):
        """Ctor.
        @param records Numpy record array with fields id,label,split. Conversion to recarray will be attempted if necessary.
        @param fileName file created by save() call
        @pre Either records or fileName must be None. id field is unique, e.g. obtained with UUID module"""
        assert records is None or fileName is None
        self.data = Struct()
        if records is not None:
            if not isinstance(records,n.ndarray):
                records = n.asarray(records,self.defRecDtype())
            self.data.records = n.unique1d(records)
            if initMaps:
                self.initMaps()
                self.check()
        elif fileName is not None:
            assert initMaps == True, "Not implemented"
            self.load(fileName)

    def __getstate__(self):
        return dict(data=self.data)

    def __setstate__(self,state):
        self.__dict__.update(state)
        self.initMaps()

    def initMaps(self):
        records = self.getRecords()
        self.idToRec = dict( (rec["id"],rec) for rec in records )
        self.labToRec = groupRecArray(records,keyField="label")
        self.splitToRec = groupRecArray(records,keyField="split")
        self.initLabToNameMap()

    @classmethod
    def defRecDtype(klass):
        """Return a default numpy dtype that can be used to construct new record arrays"""
        #idDtype comes from UUID
        return klass._defRecDtype

    @classmethod
    def makeRecords(klass,nrec,dtype=None):
        """Return a zeroed out new record array of a given size.
        The caller is responsible for filling it with values.
        @param nrec size of returned array
        @param dtype numpy dtype for array records, if None, the result of defRecDtype() will be used
        """
        if dtype is None:
            dtype = klass.defRecDtype()
        return n.zeros(nrec,dtype=dtype)

    @classmethod
    def fromIdFile(klass,featFile,label=0,split=0):
        """Create new object from IDs extracted from data file setting other fields to constant values.
        @param featFile the name of the data file. The IDs will be loaded from the accompanying id file."""
        ids = loadSeqsIdDef(featFile)
        recs = klass.makeRecords(len(ids))
        recs["id"] = ids
        recs["label"] = label
        recs["split"] = split
        return klass(records=recs)
    
    @classmethod
    def fromFields(klass,id,label=0,split=0):
        """Create new object from ID and other field arrays.
        @param id ID array
        @param label Label array (or a constant value)
        @param split Split array (or a constant value)"""
        recs = klass.makeRecords(len(id))
        recs["id"] = id
        recs["label"] = label
        recs["split"] = split
        return klass(records=recs)
    
    @classmethod
    def union(klass,idLabs):
        """Return a union of several IdLabels objects.
        @param idlabs a sequence of IdLabels objects.
        @return IdLabels object that is a union.
        @pre if there are identical ID values in input arrays, the corresponding records must be fully identical too."""
        return klass(records=n.concatenate([ x.getRecords() for x in idLabs ]))
    
    def initLabToNameMap(self):
        if not hasattr(self.data,"labNames"):
            self.data.labNames = dict( ( (lab,lab) for lab in self.labToRec ) )
        self.labToName = self.data.labNames

    def check(self):
        records = self.getRecords()
        idRecs = groupRecArray(records,keyField="id")
        nonUniqueIds = dict( ( (key,val) for (key,val) in idRecs.iteritems() if len(val) > 1 ) )
        assert len(nonUniqueIds) == 0
        
    def getRecords(self):
        return self.data.records

    def __len__(self):
        return len(self.getRecords())

    def getIdToRec(self):
        return self.idToRec

    def getIdToLab(self):
        return dict( ( (rec["id"],rec["label"]) for rec in self.getRecords() ) )

    def getLabToRec(self):
        return self.labToRec

    def getSplitToRec(self):
        return self.splitToRec

    def getSplits(self):
        """Return dict(split->IdLabel)"""
        splits = {}
        splitToRec = self.getSplitToRec()
        for key,rec in splitToRec.iteritems():
            splits[key] = self.__class__(records=rec)
        return splits

    def join(self,other,colOther="split",requireMatch=False):
        """Perform inner join with another IdLabel object or a record array.
        @param other IdLabel instance or record array with at least fields "id" and colOther
        @param colOther Data column to take from other - the rest will be taken from self
        @param requireMatch If True, will expect every ID from self to be present in other
        @return new IdLabels instance resulting from join operation"""
        if isinstance(other,IdLabels):
            idToRecOther = other.getIdToRec()
        else:
            ffarg = {}
            ffarg["id"] = other["id"]
            ffarg[colOther] = other[colOther]
            idToRecOther = self.fromFields(**ffarg).getIdToRec()
        if not requireMatch:
            recs = self.selById(set(idToRecOther),rawRec=True)
        else:
            recs = self.getRecords()
        for rec in recs:
            rec[colOther] = idToRecOther[rec["id"]][colOther]
        return self.__class__(records=recs)

    def labSet(self):
        return set(self.getLabToRec())

    def splitSet(self):
        return set(self.getSplitToRec())

    def selBySplits(self,splits):
        """Return a new IdLabels object that contains only records with 'split' value within a given list
        @param splits sequence of allowed split values
        @return new IdLabels object"""
        splits = set(splits)
        recs = self.getRecords()
        ind = n.asarray([ ind for (ind,split) in it.izip(it.count(),recs["split"]) if split in splits ],dtype=int)
        return self.__class__(records=recs[ind])

    def selByLabs(self,labs):
        """Return a new IdLabels object that contains only records within a given labels sequence.
        @todo write a single selByField() method and call that from other selByXXX()."""
        labs = set(labs)
        recs = self.getRecords()
        ind = n.asarray([ ind for (ind,lab) in it.izip(it.count(),recs["label"]) if lab in labs ],dtype=int)
        return self.__class__(records=recs[ind])

    def selById(self,ids,rawRec=False):
        """Return a new IdLabels object (or only records array) that contains only records within a given ids sequence."""
        ids = set(ids)
        recs = self.getRecords()
        ind = n.asarray([ ind for (ind,id) in it.izip(it.count(),recs["id"]) if id in ids ],dtype=int)
        records=recs[ind]
        if rawRec:
            return records
        else:
            return self.__class__(records=records)

    def selByCondition(self,condition):
        """Return a new IdLabels object that contains only records for which 'condition(record)' returns True. 
        @param condition unary predicate that will be applied to every record
        @return new IdLabels object"""
        recs = self.getRecords()
        ind = n.asarray([ ind for (ind,rec) in it.izip(it.count(),recs) if condition(rec) ],dtype=int)
        return self.__class__(records=recs[ind])
    
    def getLabToName(self):
        return self.labToName

    def getLabNames(self):
        return self.data.labNames

    def setLabNames(self,labNames):
        self.data.labNames = labNames
        self.initLabToNameMap()

    def save(self,fileName):
        dumpObj(self,fileName)

    def load(self,fileName):
        self.__dict__.update(loadObj(fileName).__dict__)

    def selDataInd(self,data):
        """Return array of row indices for data rows that are referenced here"""
        idToRec = self.getIdToRec()
        return n.asarray([ ind for (ind,id) in it.izip(it.count(),data["id"]) if id in idToRec ],dtype=int)

    def selData(self,data,setLab=True):
        """Return records from data array that are referenced here.
        @param data feature data array
        @param setLab if True, labels if returned data records will be updated from this object"""
        d = data[self.selDataInd(data)]
        if setLab:
            self.setDataLab(d)
        return d

    def setDataLab(self,data):
        """Update labels in the data feature array from this object"""
        idToRec = self.getIdToRec()
        for rec in data:
            rec["label"] = idToRec[rec["id"]]["label"]

    def balance(self,maxCount=0,targets={},rndSeed=None):
        """Return a new Idlabels object that has id counts balanced within each split across labels.
        @param maxCount default maximum count value for each (split,label) combination, @see balance() 
        at the module level for the meaning of zero and negative values.
        @param targets dict {split : { "maxCount" : value, "labTargets" : value },...} where each value is optional
        and overrides the maxCount parameter for a specific split or specific (split,label) pair.
        @return new IdLabels balanced object."""
        spRecs = self.getSplitToRec()
        spRecsBal = []
        for split in spRecs:
            spTargets = targets.get(split,{})
            spMaxCount = spTargets.get("maxCount",maxCount)
            spLabTargets = spTargets.get("labTargets",{})
            spRecsBal.append(balance(data=spRecs[split],maxCount=spMaxCount,labTargets=spLabTargets,rndSeed=rndSeed))
        return self.__class__(records=n.concatenate(spRecsBal))

    def selAcrossSplits(self):
        """Return a new IdLabels object such that each label left is present in every split.
        Rational: if we do cross-validation and label is not present in training split, but
        present in testing, we will get zero true positives for that label during testing.
        Therefore, it makes sense to do CV only when every label is present in every split.
        Alternatively we could skip during testing all samples with labels not present in
        training, but would negatively bias our specificity estimate (because corresponding
        classes would have a chance of showing false positives only."""
        splits = self.getSplits()
        labcross = set(self.getLabToRec())
        for split in splits.values():
            labcross &= set(split.getLabToRec())
        return self.selByLabs(labcross)

    def count(self,keyFields=("split","label"),format="list"):
        return countRecArray(arr=self.getRecords(),keyFields=keyFields,format=format)

    def remapIds(self,idMap):
        """Return a new IdLabels object that has new id values created from IdMap argument.
        Because IdMap represents a 1-to-N old-to-new relation, the returned object will have
        multiple copies of the same records but with different new id value."""
        oldToRec = self.getIdToRec()
        #some old ids from idMap may be absent in self, so drop them
        idMap = idMap.selByOldId(oldToRec)
        newToOld = idMap.getIdToRec()
        idmapRecs = idMap.getRecords()
        newrecs = self.makeRecords(nrec=len(idmapRecs))
        newrecs["id"] = idmapRecs["id"]
        for iRec in xrange(len(newrecs)):
            newid = newrecs[iRec]["id"]
            newrecs[iRec] = oldToRec[newToOld[newid]["oldId"]]
            newrecs[iRec] ["id"] = newid
        return self.__class__(records=newrecs)

    def renumSplits(self,start=1):
        recsO = self.getRecords()
        recsN = recsO.copy()
        splsO = recsO["split"]
        splsN = recsN["split"]
        splits = sorted(self.splitSet())
        for (i,spO) in enumerate(splits):
            splsN[splsO==spO] = i + start
        return self.__class__(records=recsN)


        
def showSvmDataCounts(data):
    return "%s records: %s" % (len(data),sorted(binCount(data['label'].astype('i4')).items()))
        
def unionIdLabels(idLabs):
    return idLabs[0].union(idLabs)

def loadIdLabelsMany(fileNames):
    return unionIdLabels([ loadObj(fileName) for fileName in fileNames ])
        

def saveIdLabelRecords(records,fileName):
    """Shortcut to saving records that can be used to build an IdLabels object as an IdLabels object w/o building internal maps.
    This avoids spending time on building IdLabels internal maps."""
    idLabs = IdLabels(records=records,initMaps=False)
    idLabs.save(fileName)


def stdScaleAndCenter(matr):
    matr = matr - matr.mean(0)
    return matr/matr.std(0)


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


def sparseDiv(x,y,alt=None):
    """return x/y where x and y are sparse features represented as (ind,val) recarays.
    @param x "ind","val" recarray
    @param y "ind","val" recarray
    @param alt "val" array with the same size as x, to use as a denominator when y index is missing
    This assumes that each feature is sorted in index order.
    @todo this is space/time inefficient for very sparse data and should be replaced by C code"""
    yd = n.zeros(y["ind"][-1]+1,dtype=y["val"].dtype)
    yd[y["ind"]] = y["val"]
    z = x.copy()
    yz = yd[z["ind"]]
    if alt is not None:
        z["val"] /= n.select([yz != 0],[yz],default=alt)
    else:
        z["val"] /= yz
    return z

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
            data = n.concatenate((data[:testRecNum/2],data[-testRecNum/2:-1]))
        kmersInp.close()
        labels = data['taxid'].astype(n.float64)
        vals = data['vals'].astype(n.float64)
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



