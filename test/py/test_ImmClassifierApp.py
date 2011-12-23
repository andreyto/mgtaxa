from MGT.ImmClassifierApp import *
from MGT.SeqDbFasta import *

jobs = []

seqDbPath1 = pjoin(options.testDataDir,"seqdb-fasta")
seqDbPath2 = pjoin(options.testDataDir,"fasta")

optTpl = Struct()
optTpl.runMode = "inproc"
optTpl.lrmUserOptions = r"'-P 0413'"

dryRun = False

cmdPref = "python $MGT_HOME/MGT/ImmClassifierApp.py"
cmdSfx = "--run-mode %s --lrm-user-options %s" %\
        (optTpl.runMode,optTpl.lrmUserOptions)

cmdLog = []

def runAndLog(cmd,help):
    cmd = "%s %s %s" % (cmdPref,cmd,cmdSfx)
    cmdLog.append((dedent(help),cmd))
    run(cmd,debug=True,dryRun=dryRun)


def buildRefSeqDbCmd():
    help="""Build the sequence DB for model training from NCBI RefSeq multi-FASTA file(s).
    This currently filters the input by excluding plasmids and taxa w/o enough total sequence."""
    cmd = "--mode make-ref-seqdb --inp-ncbi-seq '%s'" % \
        (pjoin(seqDbPath1,"*.fasta.gz"),)
    cmd += " --db-seq tmp.db-seq"
    runAndLog(cmd,help)

def trainRefCmd():
    help="Train models based on a sequence DB built by make-ref-seqdb step."
    cmd = "--mode train --db-seq tmp.db-seq --db-imm tmp.imm "+\
            "--train-max-len-samp-model 1000000 --incremental-work 1"
    runAndLog(cmd,help)

def predictAgainstRefCmd():
    help="""Make a prediction for each sequence in the --inp-seq multi-FASTA file
    against the --db-imm database of models. The output results are stored in
    --pred-out-dir. Per-sequence predictions are stored in a CSV file. Aggregated
    counts per clade at various taxonomic levels are provided in the stats sub-directory,
    along with the auto-generated graphs. The same data is provided in a SQLite file.
    Several extra options can change default locations of the individual output files.
    If you omit the --db-imm option, the program will try to use a central DB of models
    configured for this installation."""
    run("zcat %s %s > tmp.pred.fasta" % tuple([ pjoin(seqDbPath1,f) \
            for f in ("100226.fasta.gz","101510.fasta.gz") ]),\
            shell=True)
    cmd = "--mode predict --inp-seq tmp.pred.fasta --db-imm tmp.imm --pred-min-len-samp 1000"
    cmd += " --reduce-scores-early 1"
    cmd += " --pred-out-dir tmp.results"
    runAndLog(cmd,help)

def makeBenchCmd():
    help="""Create benchmark dataset based on a sequence DB built by make-ref-seqdb step.
    It also uses information from the model database built by train step."""
    cmd = "--mode make-bench  --db-seq tmp.db-seq --db-imm tmp.imm "+\
            "--db-bench tmp.db-bench --db-bench-frag tmp.bench.fna "+\
            "--db-bench-frag-size 400 --db-bench-frag-count-max 100"
    runAndLog(cmd,help)

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
    
    opt.web = True
    opt.needTerminator = True
    opt.stdout = "stdout.log"
    opt.stderr = "stderr.log"

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

#buildRefSeqDbCmd()
#trainRefCmd()
makeBenchCmd()
#predictAgainstRefCmd()
#jobs = trainRef(jobs)
#jobs = scoreRefAgainstRef(jobs)
#jobs = procScoresRefAgainstRef(jobs)

#jobs = trainCustom(jobs)
#print jobs
#jobs = scoreRefAgainstCustom(jobs)
#print jobs
#jobs = procScoresRefAgainstCustom(jobs)
#print jobs
#jobs = scoreCustomAgainstJoint(jobs)
#jobs = procScoresCustomAgainstJoint(jobs)
for x in cmdLog:
    print x[0]
    print x[1]

