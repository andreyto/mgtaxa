"""Classes for generating Makeflow inputs"""
from MGT.Common import *

class MakeflowWriter(object):
    
    def __init__(self,out,vars=None,exports=None,mode="w"):
        if isinstance(out,str):
            self.out = open(out,mode)
        #store IDs of already processed MGT jobs to avoid submitting any job twice
        self.done = set()

    def appendJob(self,targets,inputs,cmd,
            vars=None):
        w = self.out.write
        w(' '.join([str(t) for t in targets])+": "+\
                ' '.join([str(t) for t in inputs])+'\n')
        if vars:
            w("\n".join(["    {}".format(v) for v in vars])+'\n')
        w("    {}\n\n".format(cmd))

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
            A
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

