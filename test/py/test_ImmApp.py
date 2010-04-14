from MGT.ImmApp import *
from MGT.SeqDbFasta import *

seqDbPath = pjoin(options.testDataDir,"seqdb-fasta")

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

imm = ImmApp(opt=opt)
jobs = imm.run()

opt.runMode = "batchDep"
opt.mode = "score"
opt.outDir = "score"
opt.immIds = immIdToSeqIds
opt.inpSeq = pjoin(seqDbPath,"195.fasta.gz")
opt.outScoreComb = "score.comb"

imm = ImmApp(opt=opt)
imm.run(depend=jobs)


