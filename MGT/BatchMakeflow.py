"""Classes for generating Makeflow inputs"""

_makeflow_rule_tpl = """\
%(targets)s: %(inputs)s
    %(cmd)s
"""

class MakeflowWriter(object):
    
    def __init__(self,out):
        self.out = open(out,"w")
        #store IDs of already processed jobs to avoid submitting any job twice
        self.done = set()

    def _writeEntry(self,targets,inputs,cmd):
        self.out.write(_makeflow_rule_tpl % dict(
            targets=' '.join(targets),
            inputs=' '.join(inputs),
            cmd=cmd))

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
            self._writeEntry(targets=targets,inputs=inputs,cmd=cmd)
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

