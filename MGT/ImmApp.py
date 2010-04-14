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

class ImmApp(App):
    """App-derived class for building collections of IMMs/ICMs and scoring against them"""

    batchDepModes = ("score","train")

    maxSeqIdLen = UUID.maxIdLen
    maxSeqPartIdLen = UUID.maxIdLen

    immSfx = ".imm"
    scoreSfx = ".score.pkl.gz"


    def init(self):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.seqDb = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
   

    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "train":
            return self.trainMany(**kw)
        elif opt.mode == "train-one":
            return self.trainOne(**kw)
        elif opt.mode == "score":
            return self.scoreMany(**kw)
        elif opt.mode == "score-one":
            return self.scoreOne(**kw)
        elif opt.mode == "combine-scores":
            return self.combineScores(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def getSeqDb(self):
        opt = self.opt
        if self.seqDb is None:
            self.seqDb = SeqDbFasta.open(opt.seqDb) #"r"
            return self.seqDb

    def trainOne(self,**kw):
        """Train and save one IMM.
        Parameters are taken from self.opt
        @param immId Assign this ID to the IMM
        @param immSeqIds List of sequence ids from seqDb
        """
        opt = self.opt
        immId = opt.immId
        immSeqIds = opt.immSeqIds
        store = self.store
        seqDb = self.getSeqDb()
        imm = Imm(path=store.getFilePath("%s%s" % (immId,self.immSfx)))
        inp = imm.train()
        seqDb.writeFasta(ids=immSeqIds,out=inp)
        inp.close()
        imm.flush()

    def trainMany(self,**kw):
        """Train many IMMs.
        Parameters are taken from self.opt
        @param immIdToSeqIds File name that contains a dict (immId->immSeqIds)
        """
        opt = self.opt
        jobs = []
        for (immId,immSeqIds) in sorted(loadObj(opt.immIdToSeqIds).items()):
            immOpt = copy(opt)
            immOpt.mode = "train-one"
            immOpt.immId = immId
            immOpt.immSeqIds = immSeqIds
            immApp = ImmApp(opt=immOpt)
            jobs += immApp.run(**kw)
        return jobs

    def scoreOne(self,**kw):
        """Score with one IMM.
        Parameters are taken from self.opt
        @param immId Score with this IMM
        @param inpSeq Name of the input multi-FASTA file to score
        @param outScore File name for output raw scores
        """
        opt = self.opt
        immId = opt.immId
        inpFastaFile = opt.inpSeq
        outScoreFile = opt.outScore
        store = self.store
        imm = Imm(path=store.getFilePath("%s%s" % (immId,self.immSfx)))
        scores = imm.score(inp=inpFastaFile)
        dumpObj(scores,outScoreFile)
        imm.flush()

    def scoreMany(self,**kw):
        """Score with many IMMs.
        Parameters are taken from self.opt
        @param immIds List of IMM IDs to score with
        @param inpSeq Name of the input multi-FASTA file to score
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        makedir(opt.outDir)
        rmf(pjoin(opt.outDir,"*"+self.scoreSfx))
        jobs = []
        for immId in opt.immIds:
            immOpt = copy(opt)
            immOpt.mode = "score-one"
            immOpt.immId = immId
            immOpt.outScore = pjoin(opt.outDir,"%s%s" % (immId,self.scoreSfx))
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
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        scores = None
        idScores = None
        for (iImm,immId) in enumerate(opt.immIds):
            inpScoreFile = pjoin(opt.outDir,"%s%s" % (immId,self.scoreSfx))
            score = loadObj(inpScoreFile)
            if scores is None:
                scores = n.resize(score["score"],(len(score),len(opt.immIds)))
                idScores = score["id"].copy()
            else:
                assert n.all(idScores == score["id"])
            scores[:,iImm] = score["score"]
        dumpObj(ImmScores(idImms=opt.immIds,idScores=idScores,scores=scores),opt.outScoreComb)
    
    def listImmIds(self):
        return list(self.store.fileNames(pattern="*"+self.immSfx,sfxStrip=self.immSfx))

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmApp)

