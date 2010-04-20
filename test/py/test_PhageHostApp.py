from MGT.PhageHostApp import *

jobs = []

workDir = os.environ["GOSII_WORK"]

opt = Struct()
opt.runMode = "inproc"
opt.immDb = pjoin(workDir,"icm")
opt.predSeq = "samp.v/samp.v.5000/samp.fas"
opt.predOutDir = "samp.v/samp.v.5000/res"
opt.predIdLab = "samp.v/samp.v.5000/idlab.pkl"

opt.mode = "perf" #"perf" #"proc-scores"

app = PhageHostApp(opt=opt)
app.run(depend=jobs)


