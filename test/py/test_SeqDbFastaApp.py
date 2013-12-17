### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from MGT.Common import *
from MGT.SeqDbFastaApp import *

inpPath = pjoin(options.testDataDir,"seqdb-fasta","*.fasta.gz")

optTpl = Struct()
optTpl.runMode = "inproc" #"batchDep"
optTpl.lrmUserOptions = r"'-P 0413'"
#optTpl.cwd = pabs("opt_work")

#makedir(optTpl.cwd)

dryRun = False

debugger = False

cmdPref = "python %s $MGT_HOME/MGT/SeqDbFastaApp.py" % ("-m pdb" if debugger else "",)
cmdSfx = "--run-mode %s --lrm-user-options %s" %\
        (optTpl.runMode,optTpl.lrmUserOptions)

cmdLog = []

def runAndLog(cmd,help):
    cmd = "%s %s %s" % (cmdPref,cmd,cmdSfx)
    cmdLog.append((dedent(help),cmd))
    run(cmd,debug=True,dryRun=dryRun)


def splitCmd():
    rmrf("tmp.seqdb-app.split")
    help="""Test splitting of input sequence set."""
    cmd = ""
    cmd += " --mode split --inp-seq '%s'" % (inpPath,)
    cmd += " --out-seq tmp.seqdb-app.split"
    cmd += " --chunk-size 1000"
    cmd += " --out-seq-db-ids tmp.seqdb-app.split.dbids"
    cmd += " --assert-out-seq 1"
    runAndLog(cmd,help)
    cmd = ""
    cmd += " --mode split --inp-seq '%s'" % ("tmp.seqdb-app.split",)
    cmd += " --inp-seq-db-ids tmp.seqdb-app.split.dbids"
    cmd += " --out-seq tmp.seqdb-app.split.2"
    cmd += " --chunk-size 100"
    cmd += " --out-seq-db-ids tmp.seqdb-app.split.dbids.2"
    cmd += " --assert-out-seq 1"
    runAndLog(cmd,help)
    cmd = ""
    cmd += " --mode split --inp-seq '%s'" % ("tmp.seqdb-app.split.2",)
    cmd += " --inp-seq-db-ids tmp.seqdb-app.split.dbids.2"
    cmd += " --out-seq tmp.seqdb-app.split.3"
    cmd += " --chunk-size 100000000"
    cmd += " --out-seq-db-ids tmp.seqdb-app.split.dbids.3"
    cmd += " --assert-out-seq 1"
    runAndLog(cmd,help)

splitCmd()

