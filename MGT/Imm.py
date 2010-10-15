"""Wrapper around Glimmer 3 Interpolated Markov Models implementation.
"""

from MGT.Common import *

class ImmScores(object):
    """Class to represent a matrix of scores for multiple Imms (columns) for multiple samples (rows)"""
    def __init__(self,idImm,idSamp,scores,lenSamp):
        self.idImm = idImm
        self.idSamp = idSamp
        self.scores = scores
        self.lenSamp = lenSamp
    
    @classmethod
    def concatImms(klass,immScores):
        """Concatenate a sequence of ImmScores objects representing different Imms for the same samples."""
        assert len(immScores) > 0, "Need at least one object in a sequence"
        assert sameArrays((x.idSamp for x in immScores))
        idImm = n.concatenate([ x.idImm for x in immScores ])
        assert isUniqueArray(idImm), "This concatenates non-intersecting sets of Imm IDs"
        idSamp = immScores[0].idSamp.copy()
        lenSamp = immScores[0].lenSamp.copy()
        scores = n.concatenate([ x.scores for x in immScores ],1) # along rows
        return klass(idImm=idImm,idSamp=idSamp,scores=scores,lenSamp=lenSamp)
        

class Imm:
    glImmBuildExe = options.glimmer3.immBuildExe
    glImmScoreExe = options.glimmer3.immScoreExe


    def __init__(self,path):
        """Ctor.
        @param path to this IMM
        """
        self.path = path
        self.p = None
        self.cmd = None

    def __del__(self):
        self.flush()

    def train(self,inp=None):
        """Train (build) the IMM.
        @param inp Input sequences FASTA file - either an open file object, a file name, or None.
        If None, this method will return an output stream object, to which the user should write FASTA text, 
        and close it when done.
        """
        self.flush()
        cmd = shlex.split(self.glImmBuildExe+" -d 10 -w 12 -p 1 %s" % self.path)
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
        inp and out cannot both be None.
        """
        self.flush()
        print "DBG: IMM scoring with %s" % (self.path,)
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
        if out is PIPE and inp is PIPE:
            raise ValueError("Both inp and out parameters cannot be None at the same time - set one to a real file/stream")
        p = Popen(cmd, stdin=inp, stdout=out, close_fds=True)
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
            if strData is not None:
                return self.parseScores(strData.splitlines())
        else:
            return p.stdin

    def flush(self):
        """Make sure that the any subprocess trully terminated.
        This is called automatically at the start of train() and score(),
        otherwise a sequence train(); score() fails on reading the IMM header
        by IMM executable, apparently because the model file is not yet flushed.
        """
        if self.p is not None:
            self.p.wait()
            if self.p.returncode:
                raise CalledProcessError(self.cmd,self.p.returncode)

    
    @classmethod
    def parseScores(klass,inp=None):
        """Parse output from score() method.
        @param inp Can be input stream, a file name string, or a sequence of string lines
        @return record array with scoring results
        """
        closeInp = False
        if isinstance(inp,str):
            inp = openCompressed(inp,"r")
            closeInp = True
        records = []
        for line in inp:
            words = line.strip().split()
            records.append((words[0],float(words[1])))
        if closeInp:
            inp.close()
        return n.asarray(records,dtype=[("id",idDtype),("score","f4")])

