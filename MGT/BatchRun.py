"""Support for running jobs in a batch environment."""
from MGT.Util import *
from MGT.Config import *

class BatchJob(Struct):
    """Description of submitted job"""
    pass

_BatchJobTemplate = \
"""#!/bin/tcsh
#$$ -hard -P $PROJECT_CODE
#$$ -l memory=${MEM}M $LLENGTH -l arch="$ARCH"
#$$ -cwd
#$$ -r n
## Submit as 'qsub -b n -S /bin/tcsh script_name'. Apparently admins changed the default value of -b to 'y'
## and by default qstat now thinks of script_name as a binary file and does not parse it for
## embedded options (09/14/07).  Shell NEEDS to be correctly specified both at the top and (?) 
## in qstat options for the user environment to be correctly sourced.
## echo "Initial environment begin"
## printenv | sort
## echo "Initial environment end"
## pstree
source $$HOME/.cshrc
## printenv | sort
hostname
uname -a
pwd
date
top -b -n 1 | head -n 15
####
"""

class BatchSubmitter(object):
    
    @classmethod
    def defaultOptions(klass):
        opts = Struct()
        # use global options
        if hasattr(options,"batchRun"):
            options.batchRun.updateOtherMissing(opts)
        opts.setdefault('MEM',2000)
        opts.setdefault('ARCH',"lx*64")
        opts.setdefault('LENGTH',"medium")
        return opts
    
    def __init__(self,**kw):
        opts = self.defaultOptions()
        self.opts = opts
        opts.update(kw)
        length = opts.pop("LENGTH")
        if length is None:
            length = ""
        else:
            length = "-l %s" % (length,)
        opts.LLENGTH = length
        self.header = varsub(_BatchJobTemplate,**opts.asDict())
        
    def submit(self,cmd,scriptName=None,cwd=None,sleepTime=0,depend=[],dryRun=False):
        """Submit a batch job.
        @param cmd Command to run - it will be inserted verbatim into a script submitted to qsub
        @param scriptName Optional prefix for batch script name
        @param cwd Optional execution directory, otherwise the current directory will be used
        @param sleepTime sleep that many seconds after submission
        @param depend list of either job IDs or BatchJob instances on which completion the current job depends upon
        @param dryRun if True, just print what would have been done, without submitting anything
        @ret BatchJob instance for submitted job"""
        ret = None
        curdir = os.getcwd()
        try:
            if cwd:
                os.chdir(cwd)
            else:
                cwd = curdir
            if scriptName is None:
                scriptName = osNameFilter(cmd)[:10]
                if scriptName == "":
                    scriptName = "bs"
            outScr,scriptName = makeTmpFile(suffix=".qsub",prefix=scriptName+'.',dir=cwd,withTime=True)
            outScr.close()
            script = self.header + cmd + '\n'
            strToFile(script,scriptName,dryRun=dryRun)
            qsubCmd = ["qsub", "-b","n","-S","/bin/tcsh"]
            #construct dependency argument if needed
            if len(depend) > 0:
                depids = [ str(dep.jobId) if isinstance(dep,BatchJob) else str(dep) for dep in depend ]
                qsubCmd.extend(["-hold_jid",','.join(depids)])
            qsubCmd.append(scriptName)
            outp = backsticks(qsubCmd,dryRun=dryRun,dryRet="your job 0 ")
            jobId = outp.lower().split("your job")[1].strip().split()[0]
            strToFile(jobId,scriptName+".jobid",dryRun=dryRun)
            if not dryRun:
                # go easy on qsub subsystem
                sleep(sleepTime)
            ret = BatchJob(jobId=jobId,scriptName=scriptName,depend=depend)
        finally:
            os.chdir(curdir)
        return ret
        
    def nQueued(self):
        jobList = backsticks(["qstat","-u",os.environ['USER']]).strip().splitlines()
        return len(jobList)

    def submitIf(self,maxQueued,**kw):
        while self.nQueued() >= maxQueued:
            sleep(60)
        self.submit(**kw)

def runBatch(cmd,scriptName=None,cwd=None,sleepTime=0,depend=[],dryRun=False,**kw):
    if not isinstance(cmd,str):
        cmd = ' '.join(cmd)
    bs = BatchSubmitter(**kw)
    return bs.submit(cmd=cmd,scriptName=scriptName,cwd=cwd,sleepTime=sleepTime,depend=depend,dryRun=dryRun)

