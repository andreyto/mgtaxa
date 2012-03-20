### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from MGT.Common import *
from MGT.ArchiveGpgApp import *

inpPath = pjoin(options.testDataDir,"seqdb-fasta")

optTpl = Struct()
optTpl.runMode = "inproc" #"batchDep"
optTpl.lrmUserOptions = r"'-P 0413'"
#optTpl.cwd = pabs("opt_work")

#makedir(optTpl.cwd)

dryRun = False

debugger = False

cmdPref = "python %s $MGT_HOME/MGT/ArchiveGpgApp.py" % ("-m pdb" if debugger else "",)
cmdSfx = "--run-mode %s --lrm-user-options %s" %\
        (optTpl.runMode,optTpl.lrmUserOptions)

gpgArgs = "-r "+options.toolGpgKeyName

cmdLog = []

def runAndLog(cmd,help):
    cmd = "%s %s %s" % (cmdPref,cmd,cmdSfx)
    cmdLog.append((dedent(help),cmd))
    run(cmd,debug=True,dryRun=dryRun)


def archiveCmd():
    help="""Test archive creation."""
    cmd = ""
    #cmd += " --gpg-args '%s'" % (gpgArgs,)
    cmd += " --mode archive --path '%s'" % (inpPath,)
    cmd += " --archive tmp.archive.gpg"
    runAndLog(cmd,help)

def extractCmd():
    help="""Test extraction."""
    cmd = ""
    #cmd += " --gpg-args '%s'" % (gpgArgs,)
    cmd += " --mode extract --path '%s'" % ("tmp.archive.gpg.extr",)
    cmd += " --archive tmp.archive.gpg"
    cmd += " --tar-args '--strip-components 1'"
    runAndLog(cmd,help)

archiveCmd()
extractCmd()

