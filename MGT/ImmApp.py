### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application for building collections of IMMs/ICMs and scoring against them"""

from MGT.Imm import *
from MGT.Taxa import *
from MGT.Svm import *
from MGT.App import *
from MGT.DirStore import *
from MGT.SeqDbFasta import *
import UUID

from MGT.Hdf import *

class ImmScores(object):
    """Class that stores and manipulates Imm scores for a set of samples and a set of Imms.
    This is an abstract class that describes the interface.
    """

    def __init__(self,fileName,mode,**kw):
        """Constructor.
        In this default constructor, all parameters are assigned as attributes 
        to the object (by reference) and used later when the actual data 
        operations are requested.
        @param fileName name of the new object in storage
        @param mode One of w|r|c with a semantic similar to the open() builtin
        """
        self.fileName = fileName
        self.mode = mode
        for key,val in kw.items():
            setattr(self,key,val)

    def appendScore(self,idImm,score):
        """Append to the dataset a score vector for one Imm.
        @param idImm Id of the Imm
        @param score numpy.recarray(names="id","score")
        """
        pass

    def close(self):
        """Semantics is similar to a close() method on the file object (can be a noop)"""
        pass
    
    def catImms(self,fileNames,idImms=None):
        """Concatenate into self a sequence of ImmScores objects representing different Imms for the same set of samples.
        @param fileNames sequence of object names to be concatenated
        @param immIds optional sequence of all immIds that are expected to be present in the inputs;
        KeyError exception will be raised if any is missing"""
        pass

class ImmScoresHdf(ImmScores):
    """Base for implementations of ImmScores which use HDF5 as storage.
    """
    
    def setLenSamp(self,lenSamp):
        g = self.getData()
        assert n.all(g.idSamp[:] == lenSamp["id"]),"Sample IDs do not match for sample length array"
        f.createArray(g, 'lenSamp', lenSamp["len"],"Sample length")
    
    def getKind(self):
        return self.getData()._v_attrs.kind
    
    def getData(self):
        return self.hdfFile.root.immScores
    
    def open(self):
        self.hdfFile = pt.openFile(self.fileName, mode = self.mode)

    def close(self):
        """Semantics is similar to a close() method on the file object (can be a noop)"""
        if hasattr(self,"hdfFile"):
            self.hdfFile.close()

class ImmScoresReducedHdf(ImmScoresHdf):
    """Implementation of ImmScores that does score reduction and uses HDF5 as storage.
    Reduction means that at every appendScore call, only the best scores and
    associated Imm IDs are retained, in order to minimize intermediate storage.
    """
    
    def appendScore(self,idImm,score):
        """Append to the dataset a score vector for one Imm.
        @param idImm Id of the Imm
        @param score numpy.recarray(names=("id","score"))
        """
        if not self.initStorageIf(idImm,score):
            #existed already with data
            g = self.getData()
            #For PyTables datasets, the indexing operator [] is never called
            #implicitely unlike for Numpy arrays, hence [:] below.
            #Note that any [] also returns a Numpy array in memory, so beware of
            #data size.
            assert n.all(score["id"] == g.idSamp[:]),"Sample IDs do not match for scores"
            whereNewIsLarger = n.where(score["score"] > g.score[:])
            g.score[whereNewIsLarger] = score["score"][whereNewIsLarger] #Do I really need rhs selection?
            g.idImm[whereNewIsLarger] = idImm


    def initStorageIf(self,idImm,score):
        if not hasattr(self,"hdfFile"):
            self.open()
            f = self.hdfFile
            g = f.root.createGroup(f.root,"immScores")
            g._v_attrs.kind = "ImmScoresReduced"
            f.createArray(g, 'idSamp', score["id"],"Sample ID")
            idImmArr = n.zeros(len(score),dtype=idDtype)
            idImmArr[:] = idImm
            f.createArray(g, 'idImm', idImmArr, "Imm ID")
            f.createArray(g, 'score', score["score"],"Score")
            return True
        else:
            return False

    
    def catImms(self,fileNames,idImmSel=None):
        """Concatenate into self a sequence of ImmScores objects representing different Imms for the same set of samples.
        @param fileNames sequence of object names to be concatenated
        @param idImmSel optional sequence of all immIds that are expected to be present in the inputs;
        KeyError exception will be raised if any is missing"""
        assert len(fileNames) > 0,"Need a non-empty sequence of input object names"
        firstFileName = fileNames[0]
        pt.copyFile(firstFileName,self.fileName, overwrite=True)
        self.mode = "r+"
        self.open()
        f = self.hdfFile
        g = f.root.immScores
        for fileName in fileNames[1:]:
            fo = pt.openFile(fileName,mode="r")
            go = fo.root.immScores
            #For PyTables datasets, the indexing operator [] is never called
            #implicitely unlike for Numpy arrays, hence [:] below.
            #Note that any [] also returns a Numpy array in memory, so beware of
            #data size.
            assert n.all(go.idSamp[:] == g.idSamp[:]),"Sample IDs do not match for scores"
            whereNewIsLarger = n.where(go.score[:] > g.score[:])
            g.score[whereNewIsLarger] = go.score[whereNewIsLarger] #Do I really need rhs selection?
            g.idImm[whereNewIsLarger] = go.idImm[whereNewIsLarger]
            if "lenSamp" not in g and "lenSamp" in go:
                go.lenSamp._f_copy(g)
            fo.close()

class ImmScoresDenseMatrixHdf(ImmScoredHdf):
    """Implementation of ImmScores that accumulates scores in a dense matrix and uses HDF5 as storage.
    """

    def appendScore(self,idImm,score):
        """Append to the dataset a score vector for one Imm.
        @param idImm Id of the Imm
        @param score numpy.recarray(names=("id","score"))
        """
        self.initStorageIf(score)
        g = self.getData()
        #For PyTables datasets, the indexing operator [] is never called
        #implicitely unlike for Numpy arrays, hence [:] below.
        #Note that any [] also returns a Numpy array in memory, so beware of
        #data size.
        assert n.all(score["id"] == g.idSamp[:]),"Sample IDs do not match for scores"
        idImmPos = self.idImmInd[idImm]
        g.score[idImmPos,:] = score["score"]


    def initStorageIf(self,score):
        if not hasattr(self,"hdfFile"):
            self.initStorage(idSamp=samp["id"],
                    idImm=n.asarray(self.idImm,dtype=idDtype),
                    dtype=score["score"].dtype)
            return True
        else:
            return False
    
    def initStorage(self,idSamp,idImm,dtype):
        self.open()
        f = self.hdfFile
        g = f.root.createGroup(f.root,"immScores")
        g._v_attrs.kind = "ImmScoresDenseMatrix"
        f.createArray(g, 'idSamp', idSamp,"Sample ID")
        f.createArray(g, 'idImm', idImm, "Imm ID")
        self.idImmInd = dict(((x[1],x[0]) for x in enumerate(idImm)))
        #Strangely, simple Array cannot be created not from data object,
        #and so it should always fit into RAM. Hence, we use CArray here.
        #We put Imms into rows, so that adding Imms is efficient.
        #During classifcation, reduction will be done along columns (samples),
        #and we will need to do in blocks of rows (or load entire matrix at once)
        #because reading by single columns will be very I/O inefficient.
        #There are various chunkshape tricks to play for making a compromise
        #between reading and writing along different dimensions, but we are
        #not trying them yet.
        f.createCArray(g, 
                'score', 
                atom=pt.Atom.from_dtype(dtype),
                shape=(len(idImm),len(score)),
                title="Score")

    def catImms(self,fileNames,idImmSel=None):
        """Concatenate into self a sequence of ImmScores objects representing different Imms for the same set of samples.
        @param fileNames sequence of object names to be concatenated
        @param idImmSel optional sequence of all immIds that are expected to be present in the inputs;
        KeyError exception will be raised if any is missing"""
        assert len(fileNames) > 0,"Need a non-empty sequence of input object names"
        idSamp = None
        lenSamp = None
        idImm = []
        dtype = None
        for (iFile,fileName) in enumerate(fileNames):
            fo = pt.openFile(fileName,mode="r")
            go = fo.root.immScores
            if iFile == 0:
                idSamp = go.idSamp[:]
                dtype = go.score.dtype
            else:
                assert n.all(go.idSamp[:] == idSamp[:]),"Sample IDs do not match for scores"
                assert dtype == go.score.dtype,"Dtypes do not match for concatenated datasets"
            idImm.append(go.idImm[:])
            fo.close()
        idImm = n.concatenate(idImm)
        assert isUniqueArray(idImm), "This method expects non-intersecting sets of Imm IDs"

        self.mode = "w"
        self.initStorage(idSamp=idSamp,
                idImm=idImm,
                dtype=dtype)

        f = self.hdfFile
        g = f.root.immScores
        startImmPos = 0
        for (iFile,fileName) in enumerate(fileNames):
            fo = pt.openFile(fileName,mode="r")
            go = fo.root.immScores
            endImmPos = startImmPos + len(go.score)
            for pos in xrange(startImmPos,endImmPos):
                g.score[pos] = go.score[pos-startImmPos]
            startImmPos = endImmPos
            if "lenSamp" not in g and "lenSamp" in go:
                go.lenSamp._f_copy(g)
            fo.close()
            g.score.flush()


def makeImmScores(opt,*l,**kw):
    """Factory method to make a proper derivative of ImmScores class based on App options.
    """
    if opt.reduceScoresEarly:
        typeDescr = "reduced"
    elif:
        typeDescr = "dense"
    else:
        raise ValueError("Unknown value of opt.reduceScoresEarly: %s" % (opt.reduceScoresEarly,))
    if typeDescr == "reduced":
        return ImmScoresReducedHdf(*l,**kw)
    elif typeDescr == "dense":
        return ImmScoresDenseMatrixHdf(*l,**kw)
    else:
        raise ValueError("Unknown 'typeDescr' value: %s" % (typeDescr,))

class ImmStore(SampStore):
    
    immSfx = ".imm"
    
    def getImmPath(self,immId):
        return self.getFilePath("%s%s" % (immId,self.immSfx))

    def listImmIds(self,iterPaths=None):
        """List IMM IDs either from this store or from the externally provided iterable"""
        return list(self.fileNames(pattern="*"+self.immSfx,sfxStrip=self.immSfx,iterPaths=iterPaths))

def makeDefaultImmSeqIds(seqDb):
    """Return a default dict(immId -> [ seqDb-Id, ... ]) to use as opt.immIdToSeqIds in ImmApp.
    This will create a single IMM for each ID in seqDb.
    @param seqDb Either instance of SeqDbFasta or path"""

    if isinstance(seqDb,str):
        seqDb = SeqDbFasta.open(seqDb)
    return dict( [ (x,[x]) for x in seqDb.getTaxaList() ] )

class ImmApp(App):
    """App-derived class for building collections of IMMs/ICMs and scoring against them"""

    batchDepModes = ("score","train")

    maxSeqIdLen = UUID.maxIdLen
    maxSeqPartIdLen = UUID.maxIdLen

    scoreSfx = ".score"

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.setIfUndef("incrementalWork",False)
        opt.setIfUndef("immDb","imm")
        opt.setIfUndef("nImmBatches",10)
        opt.setIfUndef("trainMaxLenSampModel",10**9/2) #1G with rev-compl
        if not opt.isUndef("immIdToSeqIds"):
            opt.setIfUndef("immIds",opt.immIdToSeqIds)
        if not opt.isUndef("outDir"):
            opt.setIfUndef("outScoreComb",pjoin(opt.outDir,"combined"+klass.scoreSfx))

    def initWork(self,**kw):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.seqDb = None #will be lazy-loaded
        #self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        self.immStore = ImmStore.open(path=self.opt.immDb)
        if opt.mode == "score":
            if opt.isUndef("immIds"):
                opt.immIds = pjoin(opt.cwd,"imm-ids.pkl")
                dumpObj(self.immStore.listImmIds(),opt.immIds)

    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "train":
            return self.trainMany(**kw)
        elif opt.mode == "train-one":
            return self.trainOne(**kw)
        elif opt.mode == "score":
            return self.scoreMany(**kw)
        elif opt.mode == "score-batch":
            return self.scoreBatch(**kw)
        elif opt.mode == "combine-scores":
            return self.combineScores(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))


    def getSeqDb(self):
        opt = self.opt
        if self.seqDb is None:
            self.seqDb = SeqDbFasta.open(opt.seqDb) #"r"
        return self.seqDb


    def getImmPath(self,immId):
        return self.immStore.getImmPath(immId)

    def getScoreBatchPath(self,scoreBatchId):
        return pjoin(self.opt.outDir,"%s%s" % (scoreBatchId,self.scoreSfx))
    
    def trainOne(self,**kw):
        """Train and save one IMM.
        Parameters are taken from self.opt
        @param immId Assign this ID to the IMM
        @param immSeqIds List of sequence ids from seqDb
        """
        opt = self.opt
        immId = opt.immId
        immSeqIds = opt.immSeqIds
        seqDb = self.getSeqDb()
        immPath = self.getImmPath(immId)
        imm = Imm(path=immPath)
        try:
            inp = imm.train()
            seqDb.writeFastaBothStrands(ids=immSeqIds,out=inp,maxLen=opt.trainMaxLenSampModel)
            inp.close()
            imm.flush()
        except:
            # remove any unfinished ICM file if anything went wrong
            rmf(immPath)
            raise

    def trainMany(self,**kw):
        """Train many IMMs.
        Parameters are taken from self.opt
        @param immIdToSeqIds File name that contains a dict (immId->immSeqIds)
        """
        opt = self.opt
        #We do not generate opt.immIdToSeqIds here even when we can because
        #in that case it is difficult to make sure when we score that all
        #IMMs have been successfuly built.
        immIdToSeqIds = loadObj(opt.immIdToSeqIds)
        if not opt.incrementalWork:
            rmfMany([ self.getImmPath(immId) for immId in immIdToSeqIds.keys() ])
        jobs = []
        for (immId,immSeqIds) in sorted(immIdToSeqIds.items()):
            #if it was not incrementalWork, output was already cleaned above
            if not os.path.exists(self.getImmPath(immId)):
                immOpt = copy(opt)
                immOpt.mode = "train-one"
                immOpt.immId = immId
                immOpt.immSeqIds = immSeqIds
                immApp = ImmApp(opt=immOpt)
                jobs += immApp.run(**kw)
        #TODO: add combiner job that validates that all models have been built
        return jobs

    def scoreBatch(self,**kw):
        """Score with a batch of several IMM.
        Parameters are taken from self.opt
        @param immIds Score with these IMMs (in memory list)
        @param inpSeq Name of the input multi-FASTA file to score
        @param outDir Directory name for output score files
        """
        opt = self.opt
        immIds = opt.immIds
        inpFastaFile = opt.inpSeq
        outScoreBatchFile = self.getScoreBatchPath(opt.scoreBatchId)
        immScores = openImmScores(opt,fileName=outScoreBatchFile,mode="w",idImm=immIds)
        for immId in immIds:
            imm = Imm(path=self.getImmPath(immId))
            scores = imm.score(inp=inpFastaFile)
            immScores.appendScore(idImm=immId,scores=scores)
            imm.flush()
        immScores.close()

    def scoreMany(self,**kw):
        """Score with many IMMs.
        Parameters are taken from self.opt
        @param immIds List of IMM IDs to score with
        @param inpSeq Name of the input multi-FASTA file to score
        @param nImmBatches Number of IMM batches (determines number of batch jobs)
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        makedir(opt.outDir)
        jobs = []
        immIds = sorted(loadObj(opt.immIds))
        if opt.incrementalWork:
            #after we switched to combining within batches, it is harder to properly
            #implement the incremental scoring. In any case, it will become moot
            #after we switch to makeflow execution
            raise ValueError("Incremental mode is currently not supported for scoring")
        else:
            rmfMany([ self.getScorePath(immId) for immId in immIds ])
        rmf(opt.outScoreComb)
        immIds = n.asarray(immIds,dtype="O")
        nrnd.shuffle(immIds) # in case we doing together proks and euks
        scoreBatchIds = [] # accumulate batch ids to pass to combiner
        batches = [ x for x in enumerate(n.array_split(immIds,min(opt.nImmBatches,len(immIds)))) ]
        
        rmfMany([ self.getScoreBatchPath(scoreBatchId) for (scoreBatchId,immIdsBatch) in batches ])

        for (scoreBatchId,immIdsBatch) in batches:
            immOpt = copy(opt)
            immOpt.mode = "score-batch"
            immOpt.immIds = immIdsBatch #now this is a sequence, not a file name
            immOpt.scoreBatchId = scoreBatchId
            immApp = ImmApp(opt=immOpt)
            jobs += immApp.run(**kw)
            scoreBatchIds.append(scoreBatchId)
        coOpt = copy(opt)
        coOpt.mode = "combine-scores"
        coOpt.scoreBatchIds = scoreBatchIds
        coApp = self.factory(opt=coOpt)
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = coApp.run(**kw)
        return jobs

    def combineScores(self,**kw):
        """Combine scores as a final stage of scoreMany().
        Parameters are taken from self.opt
        @param immIds List of IMM IDs to score with
        @param inpSeq Name of the input multi-FASTA file that was scored (to pull seq lengths here)
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        scores = None
        idSamp = None
        immIds = sorted(loadObj(opt.immIds))
        immScores = openImmScores(opt,fileName=outScoreComb,mode="w")
        outScoreBatchFiles = [ self.getScoreBatchPath(scoreBatchId) \
                for scoreBatchId in opt.scoreBatchIds ]
        immScores.catImms(fileNames=outScoreBatchFiles,idImmSel=immIds)
        #@todo This will need an overhaul - we have to convert the input FASTA
        #to random access format first, and filter by length at that time, before
        #we score - otherwise we are wasting time scoring short sequences.
        lenSamp = fastaLengths(opt.inpSeq)
        immScores.setLenSamp(lenSamp)
        immScores.close()
    

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmApp)

