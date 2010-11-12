from MGT.Proj.PhHostGosApp import *

jobs = []

workDir = os.environ["GOSII_WORK"]

stage = "gos"

if stage == "ref":

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

elif stage == "gos":

    opt = Struct()
    opt.runMode = "inproc" #"inproc" #"batchDep"
    modes = ["proc-scores"] #"score-imms-gos"] #"train-imms-gos"] #"make-custom-seq"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = PhHostGosApp(opt=opt)
        jobs = app.run(depend=jobs)

else:
    raise ValueError("Unknown stage value: %s" % stage)

