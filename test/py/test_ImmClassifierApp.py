from MGT.ImmClassifierApp import *
from MGT.SeqDbFasta import *

jobs = []

seqDbPath1 = pjoin(options.testDataDir,"seqdb-fasta")
seqDbPath2 = pjoin(options.testDataDir,"fasta")

optTpl = Struct()
optTpl.runMode = "batchDep"
optTpl.PROJECT_CODE = "0413"
optTpl.LENGTH = "medium"


def trainRef(jobs):

    opt = optTpl.copy()
    opt.mode = "train"
    opt.immDb = pjoin(os.getcwd(),"imm")
    opt.seqDb = seqDbPath1

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    return jobs

def trainCustom(jobs):

    opt = optTpl.copy()
    opt.mode = "train"
    opt.inpTrainSeq = pjoin(seqDbPath2,"92830.fasta.gz")
    opt.seqDb = pjoin(os.getcwd(),"92830.seqdb")
    opt.taxaTreePkl = pjoin(os.getcwd(),"92830.tree.pkl")
    opt.immDbArchive = pjoin(os.getcwd(),"92830.immdb.tar")
    opt.trainMinLenSamp = 1

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    return jobs

def scoreRefAgainstCustom(jobs):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDbArchive = pjoin(os.getcwd(),"92830.immdb.tar")
    opt.inpSeq = pjoin(seqDbPath1,"195.fasta.gz")
    opt.outScoreComb = pjoin(os.getcwd(),"92830.combined.score.pkl.gz")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    return jobs

def scoreCustomAgainstJoint(jobs):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDb = pjoin(os.getcwd(),"imm")
    opt.immDbArchive = pjoin(os.getcwd(),"92830.immdb.tar")
    opt.inpSeq = pjoin(seqDbPath2,"92830.fasta.gz") #pjoin(seqDbPath1,"195.fasta.gz")
    opt.outScoreComb = pjoin(os.getcwd(),"92830.1.join.combined.score.pkl.gz")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    return jobs

def procScoresCustomAgainstJoint(jobs):

    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.taxaTreePkl = pjoin(os.getcwd(),"92830.tree.pkl")
    opt.outScoreComb = pjoin(os.getcwd(),"92830.1.join.combined.score.pkl.gz")
    opt.predOutDir = pjoin(os.getcwd(),"92830.1.join.results")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    return jobs

def procScoresRefAgainstCustom(jobs):

    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.taxaTreePkl = pjoin(os.getcwd(),"92830.tree.pkl")
    opt.outScoreComb = pjoin(os.getcwd(),"92830.combined.score.pkl.gz")
    opt.predOutDir = pjoin(os.getcwd(),"92830.results")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    return jobs

def scoreRefAgainstRef(jobs):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDb = pjoin(os.getcwd(),"imm")
    opt.inpSeq = pjoin(seqDbPath1,"195.fasta.gz")

    imm = ImmClassifierApp(opt=opt)
    imm.run(depend=jobs)
    return jobs

def procScoresRefAgainstRef(jobs):
    
    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.outScoreComb = pjoin(os.getcwd(),"results","combined.score.pkl.gz")

    imm = ImmClassifierApp(opt=opt)
    imm.run(depend=jobs)
    return jobs

#jobs = trainRef(jobs)
#jobs = scoreRefAgainstRef(jobs)
#jobs = procScoresRefAgainstRef(jobs)

#jobs = trainCustom(jobs)
#print jobs
#jobs = scoreRefAgainstCustom(jobs)
#print jobs
#jobs = procScoresRefAgainstCustom(jobs)
#print jobs
jobs = scoreCustomAgainstJoint(jobs)
jobs = procScoresCustomAgainstJoint(jobs)

