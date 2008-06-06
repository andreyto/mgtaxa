from types import *
import re
import string, os
from time import sleep
import time
from copy import copy, deepcopy
from cPickle import dump, load
from cStringIO import StringIO
import numpy
import numpy.random as nrnd
from tempfile import mkstemp
from textwrap import dedent

from subprocess import Popen, call, PIPE

defineCalledProcessError = False
try:
    from subprocess import CalledProcessError
except ImportError:
    defineCalledProcessError = True

if defineCalledProcessError:

    class CalledProcessError(OSError):
        def __init__(self,returncode,*l,**kw):
            OSError.__init__(self,*l,**kw)
            self.returncode = returncode

def dumpObj(obj,fileName):
    out = openCompressed(fileName,'w')
    dump(obj,out,-1)
    out.close()
    
def loadObj(fileName):
    inp = openCompressed(fileName,'rb')
    ret = load(inp)
    inp.close()
    return ret

def allChr():
    """Return a string with all characters in C local [0-255]"""
    return ''.join([chr(i) for i in range(256)])

class objectDiskCacher:

    def __init__(self,factory,fileName,recreate=False):
        self.factory = factory
        self.fileName = fileName
        self.recreate = recreate
        
    def __call__(self,*l,**kw):
        if ( not os.path.isfile(self.fileName) ) or self.recreate:
            o = self.factory(*l,**kw)
            dumpObj(o,self.fileName)
        else:
            o = loadObj(self.fileName)
        return o

def run(*popenargs, **kwargs):
    kw = {}
    kw.update(kwargs)
    dryRun = False
    if 'dryRun' in kw:
        dryRun = kw['dryRun']
        del kw['dryRun']
    if dryRun:
        print popenargs
    else:
        returncode = call(*popenargs,**kw)
        if returncode != 0:
            raise CalledProcessError(returncode=returncode)

def backsticks(*popenargs,**kwargs):
    """Similar to shell backsticks, e.g. a = `ls -1` <=> a = backsticks(['ls','-1']).
    If 'dryRun=True' is given as keyword argument, then 'dryRet' keyword must provide a value
    to return from this function."""
    kw = {}
    kw.update(kwargs)
    dryRun = False
    if 'dryRun' in kw:
        dryRun = kw['dryRun']
        del kw['dryRun']
    dryRet = None
    if 'dryRet' in kw:
        dryRet = kw['dryRet']
        del kw['dryRet']
    if dryRun:
        print popenargs
        return dryRet
    else:
        kw['stdout'] = subprocess.PIPE
        p = Popen(*popenargs, **kw)
        retout = p.communicate()[0]
        if p.returncode != 0:
            raise CalledProcessError(returncode=returncode)
        return retout


def varsub(template,*l,**kw):
    o = {}
    p = []
    for x in l:
        if isinstance(x,Struct):
            o.update(x.asDict())
        else:
            p.append(x)
    o.update(kw)
    return string.Template(template).substitute(*p,**o)

def makeTmpFile(*l,**kw):
    """Create and open a temporary file that will exist after this program closes it.
    Return a tuple (file object,file name).
    It does the same as tempfile.NamedTemporaryFile but the file is not automatically
    deleted after being closed. Because it works through calls to mkstemp and os.fdopen,
    the returned file object does not have a file name in its 'name' attribute.
    @param createParents - if True (default) - create parent directories (require 'dir' option)
    @param dir - create file in this directory
    @param mode (default 'w') - open file in this mode
    @param bufsize - open with his buffer size"""
    
    opts1 = {}
    opts1.update(kw)
    opts1.setdefault("createParents",True)
    if opts1.pop("createParents"):
        try:
            dirName = opts1["dir"]
        except KeyError:
            raise ValueError("makeTmpFile: 'dir' keyword must be used with 'createParents' keyword")
        makedir(dirName)
    l2 = []
    opts1.setdefault("mode","w")
    for k in ("mode","bufsize"):
        if opts1.has_key(k):
            l2.append(opts1[k])
            del opts1[k]
    (fd,name) = mkstemp(*l,**opts1)
    return (os.fdopen(fd,*l2),name)

def makedir(path,dryRun=False):
    run(["mkdir","-p",path],dryRun=dryRun)

#perhaps use shutil.rmtree instead?    
def rmdir(path,dryRun=False):
    run(["rm","-rf",path],dryRun=dryRun)

rmrf = rmdir

def rmf(path,dryRun=False):
    try:
        os.remove(path)
    except OSError:
        pass

def chmod(path,mode,opts='',dryRun=False):
    if isinstance(path,basestring):
        path = [path]
    else:
        path = list(path)
    run(["chmod"]+opts.split()+[mode]+path,dryRun=dryRun)

def strToFile(s,fileName,mode="w",dryRun=False):
    if not dryRun:
        out = open(fileName,mode)
        out.write(s)
        out.close()
    else:
        print "Write to file %s mode %s:" % (fileName,mode)
        print s

def openCompressed(filename,mode,**kw):
    if filename.endswith('.gz'):
        return openGzip(filename,mode,**kw)
    else:
        return open(filename,mode,2**20,**kw)

def openGzip(filename,mode,compresslevel=6):
    compresslevel = int(compresslevel)
    if mode in ("w","wb"):
        return Popen("gzip -%s > %s" % (compresslevel,filename), shell=True, env=os.environ, bufsize=2**16, stdin=PIPE, close_fds=True).stdin
    elif mode in ("r","rb"):
        return Popen(("gzip -cd %s" % (filename,)).split(),env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
    else:
        raise ValueError("'openGzip()' - Unsupported 'mode': " + mode)


def strAttributes(o,exclude=tuple(),delim=' | '):
    """Return a string with all attributes names and value for an object 'o'."""
    return delim.join([ name + ' : ' + str(val) for (name,val) in o.__dict__.items()
        if not name in exclude])

class Struct(object):
    """Class to create 'struct's on the fly.
    Example: o = Struct()
             o.i = 2
             o.x = 'ababab'
             a = Struct({'i':2,'x':'ababab'})
             b = Struct(i=2,x='ababab')
             In all three cases, the result will be the same.
             __str__ method is redefined so that printing the object of this type
             will show all attributes which are not represented as 'instance at 0xXXXXXX'.
             """

    strStyle = "p" # p - print "pretty", s - print as one string

    def __init__(self,*lw,**kw):
        for dict in lw:
            for key in dict.keys():
                setattr(self,key,dict[key])
        for key in kw.keys():
            setattr(self,key,kw[key])

    def __str__(self):
        if self.strStyle == "p":
            return self.strPretty()
        else:
            return self.strDense()

    def __getitem__(self,key):
        try:
            return getattr(self,key)
        except AttributeError:
            raise KeyError(key)

    def setdefault(self,*l):
        try:
            return getattr(self,l[0])
        except AttributeError:
            if len(l) >= 2:
                setattr(self,l[0],l[1])
                return getattr(self,l[0])
            else:
                raise KeyError(l[0])

    def update(self,other):
        if isinstance(other,Struct):
            o = other.__dict__
        else:
            o = other
        self.__dict__.update(o)
        
    def asDict(self):
        return self.__dict__

    def keys(self):
        return self.__dict__.keys()
    
    def has_key(self,key):
        return self.__dict__.has_key(key)

    def get(self,*l):
        try:
            return getattr(self,*l)
        except AttributeError:
            raise KeyError(l[0])
            
    def strDense(self):
        keys = self.keys()
        keys.sort()
        pairs = []
        for key in keys:
            obj = self.__dict__[key]
            s_obj = str(obj)
            if self.isPrintable(s_obj):
                pairs.append((key,s_obj))
        return 'Struct('+`pairs`+')'

    def __repr__(self):
        return self.__str__()

    def isPrintable(self,reprObj):
        return not re.match("^\<.+ instance at 0x[0-9a-z]+\>$",reprObj)

    def strPretty(self):
        keys = self.keys()
        keys.sort()
        s = '\n'
        for key in keys:
            obj = self.__dict__[key]
            s_obj = str(obj)
            if self.isPrintable(s_obj):
                # add to TAB to all rows of attribute's representation
                lines = s_obj.split('\n')
                s_obj = '\n\t'.join(lines)
                s = s + key + '\t=\t' + s_obj + '\n'
        return s

    def scalars(self):
        """Return dictionary mapping names of "scalar" attributes to values.
        "Scalar" attributes are non-sequence primitive types, such as Int, Float, String, None."""
        r = {}
        for key in self.keys():
            val = self.__dict__(key)
            if type(val) in (NoneType,BooleanType,IntType,LongType,FloatType,StringType):
                r[key] = val
        return r

    def copy(self):
        return copy(self)

class Options(Struct):
    
    def copy(self):
        """Deep copy semantics"""
        return deepcopy(self)

    def keys(self):
        """Will ignore all attributes that start with _"""
        return [ k for k in Struct.keys(self) if not k.startswith("_") ]

    def freeze(self):
        """Make this object read-only"""
        Struct.__setattr__(self,"_is_frozen",True)
        for name in self.keys():
            val = getattr(self,name)
            if isinstance(val,Options):
                val.freeze()

    def unfreeze(self):
        """Make this object mutable again after previous call to freeze()"""
        try:
            Struct.__delattr__(self,"_is_frozen")
        except AttributeError:
            pass
        for name in self.keys():
            val = getattr(self,name)
            if isinstance(val,Options):
                val.unfreeze()

    def __setattr__(self,name,value):
        if getattr(self,"_is_frozen",False):
            raise AttributeError(name)
        else:
            Struct.__setattr__(self,name,value)

    def __delattr__(self,name):
        if getattr(self,"_is_frozen",False):
            raise AttributeError(name)
        else:
            Struct.__delattr__(self,name,value)


class FastaReaderSink(object):
    
    def __call__(self,lineSeq):
        pass
    
class FastaReaderSinkMem(FastaReaderSink):
    
    def __init__(self):
        from cStringIO import StringIO
        self.seq = StringIO()
        
    def __call__(self,lineSeq):
        self.seq.write(lineSeq.rstrip("\n"))
        
    def sequence(self):
        return self.seq.getvalue()
    
    def close(self):
        self.seq.close()


from cStringIO import StringIO

class FastaReader(object):
    def __init__(self,infile):
        self.infile = infile
        self.freshHdr = False
        
    def records(self):
        infile = self.infile
        while True:
            if self.freshHdr:
                self.freshHdr = False
                yield self
                continue
            line = infile.readline()
            if not line:
                return
            # skip blank lines
            elif line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                yield self
    
    def header(self):
        return self.hdr

    def getNCBI_Id(self):
        """Assume that header starts with '>gi|1234567|' and return the id from second field."""
        return self.hdr.split('|',2)[1]
    
    def seqLines(self):
        infile = self.infile
        while True:
            line = infile.readline()
            if not line:
                break
            elif line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                self.freshHdr = True
                return
            yield line

    def seqChunks(self,chunkSize):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
            if seq.tell() >= chunkSize:
                yield seq.getvalue()
                seq.close()
                seq = StringIO()
        if seq.tell() > 0:
            yield seq.getvalue()
        seq.close()

    def seqArrays(self,chunkSize):
        for s in self.seqChunks(chunkSize):
            yield numpy.fromstring(s,dtype='S1')

    def sequence(self):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
        s = seq.getvalue()
        seq.close()
        return s

    def seqLen(self):
        n = 0
        for line in self.seqLines():
            n += len(line) - 1
            if not line.endswith("\n"):
                n += 1
        return n

    def close(self):
        self.infile.close()


def fastaReaderGzip(fileName):
    return FastaReader(openGzip(fileName,'r'))

def seqIterLengths(recIter):
    return numpy.fromiter((rec.seqLen() for rec in recIter),int)

def seqIterLengthsHistogram(recIter,*l,**kw):
    return numpy.histogram(seqIterLengths(recIter),*l,**kw)

class FastaRecord(object):
    def __init__(self, title, sequence):
        """'sequence' can be either a string or an integer length"""
        self.title = title
        self.sequence = sequence
        
    def __str__(self):
        hdr = '>'+self.title
        if self.hasSeq():
            return  hdr + '\n' + \
                '\n'.join([ self.sequence[i:i+70] for i in range(0,len(self.sequence),70) ])
        else:
            return hdr
    
    def hasSeq(self):
        return isinstance(self.sequence,str)
    
    def seqLen(self):
        if self.hasSeq():
            return len(self.sequence)
        else:
            return int(self.sequence)

def writeSeqByLines(out,seq,lineLen=70):
    for i in range(0,len(seq),lineLen):
        out.write(seq[i:i+lineLen])
        out.write("\n")

def readFastaRecords(infile,readSeq=True):
    saved = None
    while 1:
        # Look for the start of a record
        if saved is not None:
            line = saved
            saved = None
        else:
            line = infile.readline()
            if not line:
                return
        
        # skip blank lines
        if line.isspace():
            continue
        
        # Double-check that it's a title line
        if not line.startswith(">"):
            raise TypeError(
                "The title line must start with a '>': %r" % line)
        
        title = line.rstrip()[1:]
        
        # Read the sequence lines until the next record, a blank
        # line, or the end of file
        sequences = []
        seqLen = 0
        
        while 1:
            line = infile.readline()
            if not line or line.isspace():
                break
            if line.startswith(">"):
                # The start of the next record
                saved = line
                break
            ln = line.rstrip("\n")
            if readSeq:
                sequences.append(ln)
            else:
                seqLen += len(ln)
        
        if readSeq:
            seq = "".join(sequences)
        else:
            seq = seqLen
        yield FastaRecord(title, seq)

def splitFastaFile(inpFile,outBase,maxChunkSize):
    inpFile = open(inpFile,'r')
    inp = readFastaRecords(inpFile)
    out = None
    iChunk = 0
    chunkSize = 0
    for rec in inp:
        recSize = len(rec.sequence)
        if out is None or chunkSize + recSize > maxChunkSize:
            if out is not None:
                out.close()
            out = outBase+'_%04d'%(iChunk,)
            out = open(out,'w')
            iChunk += 1
            chunkSize = 0
        out.write('%s\n'%(rec,))
        chunkSize += recSize
    out.close()
    inpFile.close()
    return iChunk


_BatchJobTemplate = \
"""#!/bin/tcsh
#$$ -hard -P $PROJECT_CODE
#$$ -l memory=${MEM}M -l msc -l arch="$ARCH"
#$$ -cwd
## Submit as 'qsub -b n -S /bin/tcsh script_name'. Apparently admins changed the default value of -b to 'y'
## and by default qstat now thinks of script_name as a binary file and does not parse it for
## embedded options (09/14/07).  Shell NEEDS to be correctly specified both at the top and (?) 
## in qstat options for the user environment to be correctly sourced.
echo "Initial environment begin"
printenv | sort
echo "Initial environment end"
pstree
source $$HOME/.cshrc
printenv | sort
hostname
uname -a
pwd
date
top -b -n 1 | head -n 15
####
"""

class BatchSubmitter(object):
    def __init__(self,**kw):
        opts = {}
        opts.update(kw)
        opts.setdefault('MEM',2000)
        opts.setdefault('ARCH',"*amd64")
        self.header = varsub(_BatchJobTemplate,**opts)
        
    def submit(self,cmd,scriptName,runDir=None,sleepTime=1,dryRun=False):
        curdir = os.getcwd()
        try:
            if runDir:
                os.chdir(runDir)
            script = self.header + cmd + '\n'
            strToFile(script,scriptName,dryRun=dryRun)
            run(["qsub", "-b","n","-S","/bin/tcsh",scriptName],dryRun=dryRun)
            if not dryRun:
                # go easy on qsub subsystem
                sleep(sleepTime)
        finally:
            os.chdir(curdir)
        

class HistogramRdnGenerator:
    """When initialized with a histogram (as returned by numpy.histogram(...,norm=False), 
    will generate random number distributed in the same way.
    See: http://mathworld.wolfram.com/DistributionFunction.html"""
    
    def __init__(self,hist):
        self.hist = hist
        D = numpy.add.accumulate(hist[0]).astype(float)
        D = D / D[-1]
        self.D = D
        X = hist[1]
        self.X = X
        X_mid = numpy.zeros(len(X),float)
        for i in range(len(X)-1):
            X_mid[i] = (X[i] + X[i+1])/2
        X_mid[-1] = X[-1] + (X[-1]-X_mid[-2])
        D_l = numpy.roll(D,1)
        D_l[0] = 0
        self.D_l = D_l
        self.X_mid = X_mid
        self.eps = numpy.finfo(float).eps
        
    def __call__(self):
        D = self.D
        D_l = self.D_l
        X = self.X
        X_mid = self.X_mid
        y = nrnd.ranf()
        i = D.searchsorted(y)
        if D[i] - D_l[i] < self.eps:
            #random point within the current bin
            return X_mid[i] + (nrnd.ranf()-0.5)*2*(X_mid[i]-X[i])
        else:
            #linear interpolation within the bin
            return X[i] + (X_mid[i] - X[i]) * 2 * (y-D_l[i])/(D[i]-D_l[i])

    def histogram(self):
        return self.hist

def test_HistogramRdnGenerator():
    nIter = 3000
    meansSample = numpy.zeros(nIter)
    meansDistr = meansSample.copy()
    for iter in range(nIter):
        sample = nrnd.ranf(100)    
        #sample = numpy.arange(1000)
        sampleHist = numpy.histogram(sample,bins=10,normed=False)
        gen = HistogramRdnGenerator(sampleHist)
        distr = numpy.array([ gen() for i in range(len(sample)) ])
        distrHist = numpy.histogram(distr,bins=sampleHist[1],normed=False)
        meansSample[iter] = numpy.mean(sample)
        meansDistr[iter] = numpy.mean(distr)
    print distrHist, numpy.mean(distr)
    print sampleHist, numpy.mean(sample)
    print "meanDiff = ", numpy.mean(meansSample-meansDistr)


def seqLengthDistribFromSample(recIter,*l,**kw):
    lenHist = seqIterLengthsHistogram(recIter,*l,**kw)
    gen = HistogramRdnGenerator(lenHist)
    return gen
