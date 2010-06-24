from MGT.ImmClassifierApp import *
from MGT.SeqDbFasta import *

jobs = []

#seqDbPath = pjoin(options.testDataDir,"seqdb-fasta")
seqDbPath = "/usr/local/projects/GOSII/atovtchi/refseq-taxa"

opt = Struct()
opt.runMode = "batchDep"
opt.seqDb = seqDbPath

opt.mode = "train"

imm = ImmClassifierApp(opt=opt)
jobs = imm.run()
sys.exit(0)

opt.mode = "predict"
opt.inpSeq = pjoin(seqDbPath,"195.fasta.gz")

imm = ImmClassifierApp(opt=opt)
imm.run(depend=jobs)


