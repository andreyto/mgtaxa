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

    scoreSfx = ".score.pkl.gz"

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

    def getScorePath(self,immId):
        return pjoin(self.opt.outDir,"%s%s" % (immId,self.scoreSfx))
    
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
        for immId in immIds:
            outScoreFile = self.getScorePath(immId)
            imm = Imm(path=self.getImmPath(immId))
            scores = imm.score(inp=inpFastaFile)
            dumpObj(scores,outScoreFile)
            imm.flush()

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
            immIds = [ immId for immId in immIds  if not os.path.exists(self.getScorePath(immId)) ]
        else:
            rmfMany([ self.getScorePath(immId) for immId in immIds ])
        rmf(opt.outScoreComb)
        immIds = n.asarray(immIds,dtype="O")
        nrnd.shuffle(immIds) # in case we doing together proks and euks
        for immIdsBatch in n.array_split(immIds,min(opt.nImmBatches,len(immIds))):
            immOpt = copy(opt)
            immOpt.mode = "score-batch"
            immOpt.immIds = immIdsBatch #now this is a sequence, not a file name
            immApp = ImmApp(opt=immOpt)
            jobs += immApp.run(**kw)
        coOpt = copy(opt)
        coOpt.mode = "combine-scores"
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
        for (iImm,immId) in enumerate(immIds):
            inpScoreFile = pjoin(opt.outDir,"%s%s" % (immId,self.scoreSfx))
            score = loadObj(inpScoreFile)
            if scores is None:
                scores = n.resize(score["score"],(len(score),len(immIds)))
                idSamp = score["id"].copy()
            else:
                assert n.all(idSamp == score["id"])
            scores[:,iImm] = score["score"]
        lenSamp = fastaLengths(opt.inpSeq)
        assert n.all(idSamp == lenSamp["id"])
        dumpObj(ImmScores(idImm=immIds,idSamp=idSamp,scores=scores,lenSamp=lenSamp["len"]),opt.outScoreComb)
    

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmApp)

