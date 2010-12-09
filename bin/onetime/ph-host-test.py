from MGT.Proj.PhHostGosApp import *
from MGT.ImmScalingApp import *

jobs = []

topWorkDir = os.environ["GOSII_WORK"]


#stage = "build-test-data"
stage = "test"

opt = Struct()
opt.runMode = "inproc"
refname = "refseq"
#refname = "gos-bac"
opt.immDb = pjoin(topWorkDir,"icm-%s" % refname)
opt.workDir = pjoin(topWorkDir,"ph-test")
opt.cwd = opt.workDir
makedir(opt.workDir)
#opt.predSeq = pjoin(opt.workDir,"asm_combined_454_large.5K.rnd.fna")
opt.predOutDir = pjoin(opt.workDir,"pred")
opt.rndScoreComb = pjoin(topWorkDir,"icm-refseq-scale-score","combined.score.pkl.gz")
opt.nImmBatches = 200
opt.predMinLenSamp = 50 #less than shred size below
opt.shredSizeVirTest = 5000

if stage == "build-test-data":

#    for mode in ("build-db-ph","sel-db-ph-pairs","shred-vir-test"):
    for mode in ("sel-db-ph-pairs","shred-vir-test"):
        opt.mode = mode
        app = PhageHostApp(opt=opt)
        jobs = app.run(depend=jobs)

elif stage == "test":
    opt.runMode = "inproc"
    opt.predSeq = pjoin(opt.workDir,"samp.v/samp.v.5000/samp.fas")
    opt.predIdLab = pjoin(opt.workDir,"samp.v/samp.v.5000/idlab.pkl")
#    for mode in ("predict","proc-scores","perf"):
    for mode in ("proc-scores","perf"):
        opt.mode = mode
        app = PhageHostApp(opt=opt)
        jobs = app.run(depend=jobs)

elif stage == "icm-scale":
    opt.outScaleDir = pjoin(topWorkDir,"icm-%s-scale" % refname)
    opt.outScoreDir = pjoin(topWorkDir,"icm-%s-scale-score" % refname)
    
    for mode in ("score",): #generate score
        opt.mode = mode
        app = ImmScalingApp(opt=opt)
        jobs = app.run(depend=jobs)
else:
    raise ValueError("Unknown stage value: %s" % stage)

