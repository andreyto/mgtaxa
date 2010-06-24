from MGT.PhageHostApp import *

jobs = []

workDir = os.environ["GOSII_WORK"]

opt = Struct()
opt.runMode = "inproc"
opt.immDb = pjoin(workDir,"icm")
opt.predSeq = "/usr/local/projects/GOSII/shannon/Indian_Ocean_Viral/asm_combined_454_large/454LargeContigs.fna"
opt.predOutDir = "asm_combined_454_large"
opt.nImmBatches = 200
opt.predMinLenSamp = 5000

for mode in ("proc-scores",):
    opt.mode = mode #"predict" "proc-scores" #"proc-scores-phymm" #"perf" #"proc-scores"
    app = PhageHostApp(opt=opt)
    jobs = app.run(depend=jobs)


