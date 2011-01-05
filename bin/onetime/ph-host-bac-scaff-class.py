"""Classify GOS scaffolds as RefSeq clades.
A temporary hack. Works through PhageHostApp code
with rejectRanksHigher="superkingdom"."""
from MGT.Proj.PhHostGosApp import *
from MGT.ImmScalingApp import *

jobs = []

topWorkDir = os.environ["GOSII_WORK"]

topPredDir = pjoin(topWorkDir,"ph-pred")
#topPredDir = pjoin(topWorkDir,"ph-pred-random-inp")

stage = "ref"

if stage == "ref":

    opt = Struct()
    opt.runMode = "batchDep"
    refname = "refseq"
    #refname = "gos-bac"
    opt.immDb = pjoin(topWorkDir,"icm-%s" % refname)
    opt.workDir = pjoin(topWorkDir,"ph-gos-bac")
    #opt.predSeq = "/usr/local/projects/GOSII/shannon/Indian_Ocean_Viral/asm_combined_454_large/454LargeContigs.fna"
    #opt.predSeq = "/usr/local/projects/moore130/M221/M221_Newb2.5_92ID/454AllContigs.fna"
    opt.predSeq = "/usr/local/projects/GOSII/dougAnalysis/nonViral40kbPlusScaffolds.fa.gz"
    #opt.predSeq = pjoin(topWorkDir,"scaff-gos-vir","asm_combined_454_large.5K.fna")
    #opt.predSeq = pjoin(opt.workDir,"asm_combined_454_large.5K.rnd.fna")
    #opt.predOutDir = pjoin(topPredDir,"asm_combined_454_large")
    opt.predOutDir = pjoin(topPredDir,"gos-bac_refseq")
    opt.rndScoreComb = pjoin(topWorkDir,"icm-%s-scale-score" % refname,"combined.score.pkl.gz")
    opt.nImmBatches = 200
    opt.predMinLenSamp = 100000
    opt.rejectRanksHigher = "superkingdom"

    for mode in ("predict","proc-scores",):
        opt.mode = mode #"predict" "proc-scores" #"proc-scores-phymm" #"perf" #"proc-scores"
        app = PhageHostApp(opt=opt)
        jobs = app.run(depend=jobs)
    sys.exit(0)
    opt.cwd = opt.workDir
    opt.outScaleDir = pjoin(topWorkDir,"icm-%s-scale" % refname)
    opt.outScoreDir = pjoin(topWorkDir,"icm-%s-scale-score" % refname)
    
    for mode in ("score",): #generate score
        opt.mode = mode
        app = ImmScalingApp(opt=opt)
        jobs = app.run(depend=jobs)


elif stage == "gos":

    opt = Struct()
    opt.runMode = "batchDep" #"inproc" #"batchDep"
    modes = ["score-imms-gos"] #"stats-pred" "export-pred" "proc-scores" "score-imms-gos"] #"train-imms-gos"] #"make-custom-seq"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = PhHostGosApp(opt=opt)
        jobs = app.run(depend=jobs)

else:
    raise ValueError("Unknown stage value: %s" % stage)

