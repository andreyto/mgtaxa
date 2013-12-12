from MGT.ImmClassifierApp import *
from MGT.SeqDbFasta import *
from subprocess import check_call

seqDbPath1 = pjoin(options.testDataDir,"seqdb-fasta")
seqDbPath2 = pjoin(options.testDataDir,"fasta")

def run_makeflow_if(opt):
    if opt.runMode == "batchDep" and opt.batchBackend == "makeflow":
        rmf(opt.workflowFile+".makeflowlog")
        check_call([
            options.makeflow.exe,
            opt.workflowFile
            ])

optTpl = Struct()
optTpl.runMode = "inproc" #"inproc" #"batchDep"
optTpl.batchBackend = "makeflow"
optTpl.lrmUserOptions = r"'-P 0413'"
optTpl.cwd = pabs("opt_work")
optTpl.web = False
optTpl.needTerminator = True

makedir(optTpl.cwd)

dryRun = False

debugger = True

cmdPref = "python %s $MGT_HOME/MGT/ImmClassifierApp.py" % ("-m pdb" if debugger else "",)
cmdSfx = "--run-mode %s --lrm-user-options %s --batch-backend %s" %\
        (optTpl.runMode,optTpl.lrmUserOptions,optTpl.batchBackend)

cmdLog = []

def runAndLog(cmd,help):
    cmd = "%s %s %s" % (cmdPref,cmd,cmdSfx)
    cmdLog.append((dedent(help),cmd))
    run(cmd,debug=True,dryRun=dryRun)
    run_makeflow_if(optTpl)


def buildRefSeqDbCmd():
    help="""Build the sequence DB for model training from NCBI RefSeq multi-FASTA file(s).
    This currently filters the input by excluding plasmids and taxa w/o enough total sequence."""
    cmd = "--mode make-ref-seqdb --inp-train-seq '%s' --inp-train-seq-format ncbi" % \
        (pjoin(seqDbPath1,"*.fasta.gz"),)
    cmd += " --db-seq tmp.db-seq"
    runAndLog(cmd,help)

def trainRefCmd():
    help="Train models based on a sequence DB built by make-ref-seqdb step."
    cmd = "--mode train --db-seq tmp.db-seq --db-imm tmp.imm "+\
            "--train-max-len-samp-model 1000000"
    runAndLog(cmd,help)

def predictAgainstRefCmd(reduceScoresEarly,predMode):
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
    #cmd = " --mode proc-scores --out-score-comb combined.score --pred-out-stats-html tmp.krona.html --inp-seq tmp.pred.fasta --db-imm tmp.imm --pred-min-len-samp 1000"
    #runAndLog(cmd,help)
    #return
    #cmd = " --mode predict --inp-seq tmp.pred.fasta --db-imm $MGT_DATA/icm-refseq --pred-min-len-samp 1000"
    cmd = " --mode predict --inp-seq tmp.pred.fasta --db-imm tmp.imm --pred-min-len-samp 1000"
    cmd += (" --pred-mode %s" % (predMode,))
    cmd += " --score-taxids-exclude-trees 100226"
    cmd += (" --reduce-scores-early %s" % (reduceScoresEarly,))
    cmd += (" --pred-out-dir tmp.results.reduce_pol_%s_pred_mode_%s"  % (reduceScoresEarly,predMode))
    runAndLog(cmd,help)

def makeBenchCmd():
    help="""Create benchmark dataset for a given fragment length based on a sequence DB built by make-ref-seqdb step.
    It also uses information from the model database built by train step."""
    cmd = "--mode make-bench  --db-seq tmp.db-seq --db-imm tmp.imm "+\
            "--db-bench tmp.db-bench "+\
            "--db-bench-frag-len 400 --db-bench-frag-count-max 100"
    runAndLog(cmd,help)

def benchOneFragLenCmd():
    help="""Create benchmark dataset for a given fragment length based on a 
    sequence DB built by make-ref-seqdb step and evaluate the benchmark performance.
    It also uses information from the model database built by train step."""
    cmd = "--mode bench-one-frag-len  --db-seq tmp.db-seq --db-imm tmp.imm "+\
            "--bench-out-dir tmp.bench_results "+\
            "--db-bench tmp.db-bench "+\
            "--db-bench-frag-len 400 --db-bench-frag-count-max 100"
    runAndLog(cmd,help)

def benchCmd():
    help="""Create benchmark dataset for a given list of fragment lengths based on a 
    sequence DB built by make-ref-seqdb step and evaluate the benchmark performance.
    It also uses information from the model database built by train step."""
    cmd = "--mode bench  --db-seq tmp.db-seq --db-imm tmp.imm "+\
            "--bench-out-dir tmp.bench_results "+\
            "--db-bench tmp.db-bench "+\
            "--bench-frag-len-list 50,400 --db-bench-frag-count-max 100"
    runAndLog(cmd,help)

def makeSeqDbRef(jobs):

    opt = optTpl.copy()
    opt.mode = "make-ref-seqdb"
    opt.inpTrainSeq = pjoin(seqDbPath1,"*.fasta.gz")
    opt.inpTrainSeqFormat = "ncbi"

    opt.immDb = [pjoin(opt.cwd,"imm")]
    opt.seqDb = pjoin(opt.cwd,"seqdb")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def trainRef(jobs):

    opt = optTpl.copy()
    opt.mode = "train"
    opt.immDb = [pjoin(opt.cwd,"imm")]
    opt.seqDb = pjoin(opt.cwd,"seqdb")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def makeSeqDbCustom(jobs):

    opt = optTpl.copy()
    opt.mode = "make-ref-seqdb"
    opt.inpTrainSeq = pjoin(seqDbPath2,"generic.mod.train.fasta.gz")
    opt.inpTrainModelDescr = pjoin(seqDbPath2,"generic.mod.train.json")
    opt.inpTrainSeqFormat = "generic"
    opt.seqDb = pjoin(opt.cwd,"92830.seqdb")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def trainCustom(jobs):

    opt = optTpl.copy()
    opt.mode = "train"
    opt.seqDb = pjoin(opt.cwd,"92830.seqdb")
    opt.immDb = [pjoin(opt.cwd,"92830.immdb")]
    opt.trainMinLenSamp = 1
    
    opt.stdout = "stdout.log"
    opt.stderr = "stderr.log"

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def trainCustomWithParent(jobs):

    opt = optTpl.copy()
    opt.mode = "train"
    opt.inpTrainSeq = pjoin(seqDbPath2,"custom_with_parent.fasta.gz")
    opt.seqDb = pjoin(opt.cwd,"custom_with_parent.seqdb")
    opt.taxaTreePkl = pjoin(opt.cwd,"custom_with_parent.tree.pkl")
    #opt.immDbArchive = [pjoin(opt.cwd,"custom_with_parent.immdb.tar")]
    opt.immDb = [pjoin(opt.cwd,"custom_with_parent.immdb")]
    opt.trainMinLenSamp = 1
    
    opt.stdout = "stdout.log"
    opt.stderr = "stderr.log"

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def scoreRefAgainstCustom(jobs,inpIsSeqDb=False):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDb = [pjoin(opt.cwd,"92830.immdb")]
    if inpIsSeqDb:
        opt.inpSeq = pjoin(opt.cwd,"92830.seqdb")
    else:
        opt.inpSeq = pjoin(seqDbPath1,"195.fasta.gz")
    opt.outScoreComb = pjoin(opt.cwd,"92830.combined.score")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def scoreCustomAgainstJoint(jobs):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDb = [pjoin(opt.cwd,"imm"), pjoin(opt.cwd,"92830.immdb")]
    opt.inpSeq = pjoin(seqDbPath2,"92830.fasta.gz") #pjoin(seqDbPath1,"195.fasta.gz")
    opt.outScoreComb = pjoin(opt.cwd,"92830.1.join.combined.score")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def scoreCustomWithParentAgainstJoint(jobs):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDb = [pjoin(opt.cwd,"imm"),pjoin(opt.cwd,"custom_with_parent.immdb")]
    opt.inpSeq = pjoin(seqDbPath2,"custom_with_parent.fasta.gz")
    opt.outScoreComb = pjoin(opt.cwd,"custom_with_parent.join.combined.score")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def procScoresCustomAgainstJoint(jobs):

    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.outScoreComb = pjoin(opt.cwd,"92830.1.join.combined.score")
    opt.predOutDir = pjoin(opt.cwd,"92830.1.join.results")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def procScoresCustomWithParentAgainstJoint(jobs):

    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.taxaTreePkl = pjoin(opt.cwd,"custom_with_parent.tree.pkl")
    opt.outScoreComb = pjoin(opt.cwd,"custom_with_parent.join.combined.score")
    opt.predOutDir = pjoin(opt.cwd,"custom_with_parent.join.results")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def procScoresRefAgainstCustom(jobs,inpIsSeqDb=False):

    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.outScoreComb = pjoin(opt.cwd,"92830.combined.score")
    opt.predOutDir = pjoin(opt.cwd,"92830.results")
    if inpIsSeqDb:
        opt.sampAttrib = None
    else:
        opt.sampAttrib = pjoin(seqDbPath1,"195.immClassifier.attrib.csv")

    ImmClassifierApp.fillWithDefaultOptions(opt)

    print opt

    imm = ImmClassifierApp(opt=opt)
    jobs = imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def scoreRefAgainstRef(jobs):

    opt = optTpl.copy()
    opt.mode = "score"
    opt.immDb = [pjoin(opt.cwd,"imm")]
    opt.inpSeq = pjoin(seqDbPath1,"195.fasta.gz")
    opt.outScoreComb = pjoin(opt.cwd,"results","combined.score")

    imm = ImmClassifierApp(opt=opt)
    imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def procScoresRefAgainstRef(jobs):
    
    opt = optTpl.copy()
    opt.mode = "proc-scores"
    opt.outScoreComb = pjoin(opt.cwd,"results","combined.score")
    opt.predOutDir = pjoin(opt.cwd,"results")

    imm = ImmClassifierApp(opt=opt)
    imm.run(depend=jobs)
    run_makeflow_if(opt)
    return jobs

def runAllTests():
    jobs = []
    buildRefSeqDbCmd()
    trainRefCmd()
    makeBenchCmd()
    benchOneFragLenCmd()
    benchCmd()
    predictAgainstRefCmd(reduceScoresEarly=0,predMode="taxa")
    predictAgainstRefCmd(reduceScoresEarly=1,predMode="taxa")
    
    jobs = makeSeqDbRef(jobs)
    jobs = trainRef(jobs)
    jobs = scoreRefAgainstRef(jobs)
    jobs = procScoresRefAgainstRef(jobs)

    jobs = trainCustom(jobs)
    print jobs
    inpIsSeqDb=True
    jobs = scoreRefAgainstCustom(jobs,inpIsSeqDb=inpIsSeqDb)
    print jobs
    jobs = procScoresRefAgainstCustom(jobs,inpIsSeqDb=inpIsSeqDb)
    print jobs
    jobs = scoreCustomAgainstJoint(jobs)
    jobs = procScoresCustomAgainstJoint(jobs)

    jobs = trainCustomWithParent(jobs)
    print jobs
    jobs = scoreCustomWithParentAgainstJoint(jobs)
    jobs = procScoresCustomWithParentAgainstJoint(jobs)

def runSomeTests():
    jobs = []
    if False:
        buildRefSeqDbCmd()
        trainRefCmd()
    if True:
        #makeBenchCmd()
        benchOneFragLenCmd()
        #benchCmd()
    if False:
        #predictAgainstRefCmd(reduceScoresEarly=0,predMode="host")
        predictAgainstRefCmd(reduceScoresEarly=1,predMode="taxa")
        
        jobs = makeSeqDbRef(jobs)
        jobs = trainRef(jobs)
        jobs = scoreRefAgainstRef(jobs)
        jobs = procScoresRefAgainstRef(jobs)
        jobs = makeSeqDbCustom(jobs)
        jobs = trainCustom(jobs)
        #print jobs
        inpIsSeqDb=True
        jobs = scoreRefAgainstCustom(jobs,inpIsSeqDb=inpIsSeqDb)
        #print jobs
        jobs = procScoresRefAgainstCustom(jobs,inpIsSeqDb=inpIsSeqDb)
        #print jobs
        jobs = scoreCustomAgainstJoint(jobs)
        jobs = procScoresCustomAgainstJoint(jobs)

        #jobs = trainCustomWithParent(jobs)
        #print jobs
        #jobs = scoreCustomWithParentAgainstJoint(jobs)
        #jobs = procScoresCustomWithParentAgainstJoint(jobs)

runSomeTests()

for x in cmdLog:
    print x[0]
    print x[1]

