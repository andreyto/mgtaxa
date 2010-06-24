from MGT.PhageHostApp import *

jobs = []

workDir = os.environ["GOSII_WORK"]

opt = Struct()
opt.runMode = "inproc"
opt.immDb = pjoin(workDir,"icm")
opt.predSeq = "samp.v/samp.v.5000/samp.fas"
opt.predOutDir = "samp.v/samp.v.5000/res"
#opt.predOutTaxa = pjoin(opt.predOutDir,"pred-taxa.phymm")
opt.predIdLab = "samp.v/samp.v.5000/idlab.pkl"
opt.outPhymm = "/usr/local/projects/GOSII/atovtchi/phymm/results.01.phymm____ph_samp_v_samp_v_5000_samp_fas.txt"

for mode in ("proc-scores","perf"):
    opt.mode = mode #"proc-scores" #"proc-scores-phymm" #"perf" #"proc-scores"
    app = PhageHostApp(opt=opt)
    jobs = app.run(depend=jobs)


