### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application for building collections of IMMs/ICMs"""

from MGT.Imm import *
from MGT.Taxa import *
from MGT.Svm import *
from MGT.App import *
from MGT.DirStore import *
import UUID
from MGT.SeqImportApp import *

class ImmApp(App):
    """App-derived class for building collections of IMMs/ICMs"""

    #batchDepModes = ("blastcr","pilercr")

    maxSeqIdLen = UUID.maxIdLen
    maxSeqPartIdLen = UUID.maxIdLen

    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        optChoicesMode = ("build-db-ph","sel-db-ph-pairs","shred-vir-test","predict")
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=optChoicesMode,
            dest="mode",default="predict",help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            make_option(None, "--db-gb-inp",
            action="store", type="string",dest="dbGbInp",help="Input viral GenBank file for phage-host DB construction"),
            make_option(None, "--db-ph",
            action="store", type="string",default="ph",dest="dbPh",help="DB of known phage-host pairs created and used here"),
            make_option(None, "--db-seq",
            action="append", type="string",dest="dbSeqInp",help="Input DB (e.g. RefSeq) FASTA sequence files for both microbes "+\
                    "and viruses referenced by phage-host pairs DB. It can be a subset of superset of p-h DB."),
            make_option(None, "--db-ph-pairs",
            action="store", type="string",default="ph-pairs",dest="dbPhPairs",help="Stratified phage-host pairs from DB phage-host"),
            make_option(None, "--db-ph-pairs-seq-ids",
            action="store", type="string",default="ph-pairs-seq-ids",dest="dbPhPairsSeqIds",help="Sequence IDs for db-ph-pairs - "+\
                    "intermediate file used to pull all sequences from the original FASTA files into db-ph-pairs-seq"),
            make_option(None, "--db-ph-pairs-seq",
            action="store", type="string",default="ph-pairs-seq",dest="dbPhPairsSeq",help="Sequences for db-ph-pairs in our "+\
                    "internal format"),
            make_option(None, "--shred-size-vir-test",
            action="store", type="int",default=400,dest="shredSizeVirTest",help="Shred DB viral samples into this size for testing"),
        ]
        return Struct(usage = "Build phage-host association database, classifiers, and perform training, validation and de novo classification.\n"+\
                "%prog [options]",option_list=option_list)

    
    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        refSeqDir = globals()["options"].refSeqDataDir
        opt.setIfUndef("dbGbInp",pjoin(refSeqDir,"viral.genomic.gbff.gz"))
        opt.setIfUndef("dbSeqInp",[ pjoin(refSeqDir,div+".genomic.fna.gz") for div in ("microbial","viral")])
   

    def init(self):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        self.tmpDir = self.store.getFilePath("tmp")
        makedir(self.tmpDir)
   

    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "build-db-ph":
            return self.buildDbPhageHost()
        elif opt.mode == "sel-db-ph-pairs":
            return self.selectPhageHostPairs()
        elif opt.mode == "shred-vir-test":
            return self.shredVirTest()
        elif opt.mode == "predict":
            return self.predict()
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

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
        imm = Imm(path=store.getFilePath("%s.imm" % immId))
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
        for (immId,immSeqIds) in sorted(opt.immIdToSeqIds.items()):
            immOpt = copy(opt)
            immOpt.mode = "train-one"
            immOpt.immId = immId
            immOpt.immSeqIds = immSeqIds
            immApp = ImmApp(opt=immOpt)
            jobs.append(immApp.run(**kw))
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
        imm = Imm(path=store.getFilePath("%s.imm" % immId))
        scores = imm.score(inp=inpFastaFile)
        dumpObj(scores,outScoreFile)
        imm.flush()

    def scoreMany(self,**kw):
        """Score with many IMMs.
        Parameters are taken from self.opt
        @param immIds List of IMM IDs to score with
        @param inpSeq Name of the input multi-FASTA file to score
        @param outDir Directory name for output score files
        """
        opt = self.opt
        makedir(opt.outDir)
        rmf(pjoin(opt.outDir,"*.score.gz"))
        jobs = []
        for immId in opt.immIds:
            immOpt = copy(opt)
            immOpt.mode = "score-one"
            immOpt.immId = immId
            immOpt.outScore = pjoin(opt.outDir,"%s.score.gz" % immId)
            immApp = ImmApp(opt=immOpt)
            jobs.append(immApp.run(**kw))
        return jobs

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmApp)

