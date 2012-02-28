from MGT.ImmApp import *
from MGT.SeqDbFasta import *

seqDbPath = pjoin(options.testDataDir,"seqdb-fasta")

jobs = []

opt = Struct()
opt.runMode = "inproc" #"batchDep"
opt.seqDb = seqDbPath
opt.immDb = pabs("test.immdb")
seqDb = SeqDbFasta(opt.seqDb)
ids = seqDb.getIdList()
for id in ids:
    seqDb.finById(id=id)
immIdToSeqIds = dict(((id,[id]) for id in ids))
immIdToSeqIdsFile = pabs("test.immapp.seqids.pkl")
dumpObj(immIdToSeqIds,immIdToSeqIdsFile)
opt.immIdToSeqIds = immIdToSeqIdsFile

opt.mode = "train"

ImmApp.fillWithDefaultOptions(opt)

#imm = ImmApp(opt=opt)
#jobs = imm.run()

for (reduceScoresEarly,cwd) in ((1,"imm.test.red_early_1"),(0,"imm.test.red_early_0")):

    immIds = ImmStore.open(path=opt.immDb,mode="r").listImmIdsWithIdScoreIdent()
    immIdsFile = pabs("test.immapp.immids.pkl")
    dumpObj(immIds,immIdsFile)
    opt.cwd = pabs(cwd)
    opt.mode = "score"
    opt.reduceScoresEarly = reduceScoresEarly
    opt.outDir = pabs(pjoin(opt.cwd,"score"))
    opt.immIds = immIdsFile
    opt.nImmBatches = 3
    opt.inpSeq = pjoin(seqDbPath,"195.fasta.gz")
    opt.outScoreComb = pabs(pjoin(opt.cwd,"score.comb"))

    imm = ImmApp(opt=opt)
    imm.run(depend=jobs)


