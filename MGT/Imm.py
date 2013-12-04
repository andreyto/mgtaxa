"""Wrapper around Glimmer 3 Interpolated Markov Models implementation.
"""

from MGT.Common import *

import tempfile

class Imm:
    glImmBuildExe = options.glimmer3.immBuildExe
    glImmScoreExe = options.glimmer3.immScoreExe

    workModelSfx = ".work"

    ## numpy.dtype for Imm score
    scoreDtype = "f4"

    def __init__(self,path):
        """Ctor.
        @param path to this IMM
        """
        self.path = path
        self.p = None
        self.cmd = None
        self.postFlush = None

    def __del__(self):
        self.flush()

    def train(self,inp=None):
        """Train (build) the IMM.
        @param inp Input sequences FASTA file - either an open file object, a file name, or None.
        If None, this method will return an output stream object, to which the user should write FASTA text, 
        and close it when done.
        """
        self.flush()
        workPath = self.path + self.workModelSfx
        cmd = shlex.split(self.glImmBuildExe+" -d 10 -w 12 -p 1 %s" % (workPath,))
        closeInp = False
        if inp is None:
            inp = PIPE
        elif isinstance(inp,str):
            inp = openCompressed(inp,"r")
            closeInp = True
        else:
            closeInp = False
        p = Popen(cmd, stdin=inp, close_fds=True)
        self.p = p
        self.cmd = cmd
        def _postFlushAction():
            if os.path.isfile(workPath):
                if os.path.getsize(workPath) == 0:
                    raise ValueError("Model work file has size zero after training: %s" % (workPath,))
                os.rename(workPath,self.path)
        self.postFlush = _postFlushAction
        if closeInp:
            inp.close()
        if inp is not PIPE:
            p.communicate()
            if self.p.returncode:
                raise CalledProcessError(self.cmd,self.p.returncode)
        else:
            return p.stdin

    
    def score(self,inp=None,out=None):
        """Score input sequences with an IMM.
        @param inp Input sequences FASTA file, either an open file object, a file name, or None.
        @param out Output scores file - either an open file object, a file name, or None.
        @return If inp is None, return output stream to which the user should write FASTA text and 
        close it when done; If out is None, return array of output records.
        If output stream is returned, call parseScores(inp=None) afterwards.
        """
        self.flush()
        cmd = [self.glImmScoreExe,"-N",self.path]
        closeInp = False
        closeOut = False
        if inp is None:
            inp = PIPE
        elif isinstance(inp,str):
            inp = openCompressed(inp,"r")
            closeInp = True
        else:
            closeInp = False
        if out is None:
            out = PIPE
        elif isinstance(out,str):
            out = openCompressed(out,"w")
            closeOut = True
        else:
            closeOut = False
        stderr = open(os.devnull,"w") #incompat with close_fds on Windows
        #There are reports that PIPE can fail or block on large outputs,
        #so we replace it here with a temp file in the current working
        #directory
        scoreTmpName = None
        if out is PIPE:
            out,scoreTmpName = tempfile.mkstemp(suffix=".tmp", prefix="imm.score.", dir=os.getcwd())
        self.scoreTmpName = scoreTmpName
        p = Popen(cmd, stdin=inp, stdout=out, stderr=stderr)
        stderr.close()
        self.p = p
        self.cmd = cmd
        if closeInp:
            inp.close()
        if closeOut:
            out.close()
        if inp is not PIPE:
            strData=p.communicate()[0]
            if self.p.returncode:
                raise CalledProcessError(self.cmd,self.p.returncode)
            if scoreTmpName is not None:
                out = os.fdopen(out,"r")
                out.seek(0)
                ret = self.parseScores(out)
                #os.close(out)
                out.close()
                os.remove(scoreTmpName)
                self.scoreTmpName = None
                return ret
        else:
            if scoreTmpName is not None:
                os.close(out)
            return p.stdin

    def flush(self):
        """Make sure that any subprocess trully terminated.
        This is called automatically at the start of train() and score(),
        otherwise a sequence train(); score() fails on reading the IMM header
        by IMM executable, apparently because the model file is not yet flushed.
        """
        if self.p is not None:
            self.p.wait()
            if self.p.returncode:
                raise CalledProcessError(self.cmd,self.p.returncode)
        if self.postFlush is not None:
            self.postFlush()

    
    def parseScores(self,inp=None,removeInp=True):
        """Parse output from score() method.
        @param inp Can be input stream, a file name string, or a sequence of string lines
        @param if True and inp is file name, remove inp after reading
        @return record array with scoring results
        """
        closeInp = False
        if inp is None:
            inp = self.scoreTmpName
            inpIsScoreTmpName = True
        else:
            inpIsScoreTmpName = False
        if is_string(inp):
            inpName = inp
            inp = openCompressed(inp,"r")
            closeInp = True
        records = []
        for line in inp:
            words = line.strip().split()
            records.append((words[0],float(words[1])))
        if closeInp:
            inp.close()
            if removeInp:
                os.remove(inpName)
                if inpIsScoreTmpName:
                    self.scoreTmpName = None
        return n.asarray(records,dtype=[("id",idDtype),("score",self.scoreDtype)])

