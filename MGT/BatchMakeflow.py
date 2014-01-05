"""Classes for generating Makeflow inputs"""
from MGT.Common import *

class MakeflowWriter(object):

    _tab = " "*4
    
    def __init__(self,out,vars=None,exports=None,mode="w"):
        if is_string(out):
            self.out = open(out,mode)
            #do not repeat default exports on append -
            #this breaks the Makeflow
            if "a" in mode and exports is None:
                exports = []
        elif exports is None:
            #do not add default exports when existing
            #stream is supplied. The caller should use
            #appendInitExports() when needed.
            exports = []

        self.appendInitExports(exports=exports)
        self._write_vars(vars,0)
        #store IDs of already processed MGT jobs to avoid submitting any job twice
        self.done = set()
    
    @classmethod
    def getStandardExports(klass):
        return [
                 "BATCH_OPTIONS",
                 "MAKEFLOW_BATCH_QUEUE_TYPE",
                 "MAKEFLOW_MAX_REMOTE_JOBS"
            ]

    def appendInitExports(self,exports):
        if exports is None:
            exports = self.getStandardExports()
        w = self.out.write
        for export in exports:
            w("export {}\n".format(export))

    def _write_vars(self,vars,tab_level=0):
        if vars is None:
            vars = []
        w = self.out.write
        for var in vars:
            w(self._tab*tab_level+"{}\n".format(var))
    
    def appendJob(self,cmd,targets=[],inputs=[],
            vars=None):
        w = self.out.write
        w(' '.join([str(t) for t in targets])+": "+\
                ' '.join([str(t) for t in inputs])+'\n')
        self._write_vars(vars,1)
        w(self._tab+"{}\n\n".format(cmd))

    def appendMakeflow(self,flow,targets=[],inputs=[],
            vars=None):
        """Append a task that itself is a Makeflow.
        This will mark the sub-makeflow as local unless
        it is already marked otherwise in the 'vars'"""
        if flow not in inputs:
            inputs = inputs + [flow]
        if vars is None:
            vars = []
        for var in vars:
            if var.strip().startswith("@BATCH_LOCAL"):
                break
        else:
            vars = vars + ["@BATCH_LOCAL=1"]
        self.appendJob(targets=targets,inputs=inputs,
                cmd="MAKEFLOW {}".format(flow),
                vars=vars)

    def appendMgtJob(self,job):
        """Append to current makeflow one MGTAXA job object.
        No recursion to append dependencies is performed.
        @param job BatchJob object, having at least these
        attributes: BatchJob(jobId,scriptName=scriptName,cwd=cwd,outputs=(flagOk,),depend=depend)
        where jobId is globally unique.
        """
        jobId = job.jobId
        if jobId not in self.done:
            inputs = []
            for dep in job.depend:
                inputs += dep.outputs
            targets = job.outputs
            cmd = "bash " + job.scriptName
            self.appendJob(targets=targets,inputs=inputs,cmd=cmd)
            self.done.add(job.jobId)

        
    def appendMgtJobs(self,jobs):
        """Append to current makeflow from a graph of MGTAXA job objects.
        This will recurse first to append dependencies.
        @param jobs a list of BatchJob objects, @see appendMgtJob for requirements.
        """
        for job in jobs:
            self.appendMgtJobs(job.depend)
            self.appendMgtJob(job)

    def close(self):
        if self.out:
            self.out.close()
            self.out = None


def writeMakeflowRunScript(makeflow,workflow,env,vars,args,out,mode="w"):
    """Write a shell script that will run this makeflow.
    @param makeflow Makeflow executable path
    @param workflow workflow file path
    @param env file to source as shell environment
    @param vars list of environment variable assignments "VAR=VAL"
    @param args string with all arguments to makeflow executable
    @param out file path or file object for writing the script into
    @param mode to open out if out if a file path"""
    out_close = False
    if is_string(out):
        out = open(out,mode)
        out_close = True
    w = out.write
    w("#!/bin/bash\n")
    w(". {}\n".format(env))
    for var in vars:
        w("export {}\n".format(var))
    w("{makeflow} {args} {workflow}".\
            format(
                makeflow = makeflow,
                args = args,
                workflow = workflow))
    if out_close:
        out.close()

