from MGT.ImmApp import *
from MGT.SeqDbFasta import *

seqDbPath = pjoin(options.testDataDir,"seqdb-fasta")

jobs = []

opt = Struct()
opt.runMode = "batchDep"
opt.seqDb = seqDbPath
seqDb = SeqDbFasta(opt.seqDb)
ids = seqDb.getIdList()
immIdToSeqIds = dict(((id,[id]) for id in ids))
immIdsFile = "test.immapp.seqids.pkl"
dumpObj(immIdToSeqIds,immIdsFile)
opt.immIdToSeqIds = immIdsFile

opt.mode = "train"

#imm = ImmApp(opt=opt)
#jobs = imm.run()

opt.mode = "score"
opt.outDir = "score"
opt.immIds = immIdsFile
opt.nImmBatches = 3
opt.inpSeq = pjoin(seqDbPath,"195.fasta.gz")
opt.outScoreComb = "score.comb"

imm = ImmApp(opt=opt)
imm.run(depend=jobs)


