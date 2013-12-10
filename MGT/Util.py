### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from __future__ import with_statement  # only for python 2.5
from contextlib import contextmanager, closing
from types import *
import re
import string, os, sys
pjoin = os.path.join
pabs = os.path.abspath
from time import sleep
import time
from copy import copy, deepcopy
from cPickle import dump, load, Pickler, Unpickler
from cStringIO import StringIO
import json
import tempfile

import numpy

n = numpy
nma = numpy.ma
import numpy.random as nrnd

np = numpy
npma = numpy.ma
import numpy.random as nprnd

import random

from tempfile import mkstemp
from textwrap import dedent
from collections import defaultdict as defdict
import operator
import itertools as it
import datetime
import shutil
import contextlib
import bz2
import pdb

from MGT.Options import Struct
from MGT.Config import options, Options

from MGT.Run import *

from MGT.Bits.Unique import unique, uniquePick

@contextlib.contextmanager
def chdir(dirname=None):
    curdir = os.getcwd()
    try:
        if dirname is not None:
            os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)


class Timer(object):
    """Very simple timer object"""

    def __init__(self):
        ##time.clock() seems to be broken on SuSe 10 x86_64 Python 2.4
        ##- it always returns the same value
        self.start = time.time()

    def __call__(self):
        finish = time.time()
        lap = finish-self.start
        self.start = finish
        return lap

    def msg(self,s):
        return "%s. Elapsed: %s" % (s,self())

def dumpObj(obj,fileName,**kw):
    #print "DEBUG: dumpObj({},{},{})".format(obj,fileName,kw)
    with closing(openCompressed(fileName,'wb',**kw)) as out:
        dump(obj,out,-1)
    
def loadObj(fileName,**kw):
    with closing(openCompressed(fileName,'rb',**kw)) as inp:
        return load(inp)

def dumpNumpy(obj,fileName,**kw):
    #numpy.save does not like non-file streams
    #out = openCompressed(fileName,'w',**kw)
    with closing(open(fileName,'wb',**kw)) as out:
        np.save(out,obj)
    
def loadNumpy(fileName,**kw):
    #inp = openCompressed(fileName,'rb',**kw)
    with closing(open(fileName,'rb',**kw)) as inp:
        return np.load(inp)

def PickleReader(inp,**kw):
    if isinstance(inp,str):
        closeInp = True
        inp = openCompressed(inp,'rb',**kw)
    else:
        closeInp = False
    pkl = Unpickler(inp)
    try:
        while True:
            try:
                x = pkl.load()
            except EOFError:
                break
            yield x
    finally:
        if closeInp:
            inp.close()

class PickleWriter(object):

    def __init__(self,out,**kw):
        self.out = out
        if isinstance(out,str):
            self.closeOut = True
            self.out = openCompressed(out,'w',**kw)
        else:
            self.closeOut = False
        self.pkl = Pickler(self.out,-1)

    def write(self,x):
        self.pkl.dump(x)
    
    def close(self):
        self.out.close()

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

def currTimeAsFileName():
    """Return current datetime formatted as suitable for use as a file name"""
    return datetime.datetime.today().strftime("%y-%m-%d_%H-%M")

def makeTmpFile(*l,**kw):
    """Create and open a temporary file that will exist after this program closes it.
    @return a tuple (file object,file name).
    It does the same as tempfile.NamedTemporaryFile but the file is not automatically
    deleted after being closed. Because it works through calls to mkstemp and os.fdopen,
    the returned file object does not have a file name in its 'name' attribute.
    @param prefix If provided, the file will start with that prefix (inside the dir directory)
    @param suffix If provided, the file will have that suffix
    @param createParents - if True (default) - create parent directories (require 'dir' option)
    @param dir - create file in this directory
    @param mode (default 'w') - open file in this mode
    @param bufsize - open with this buffer size"""
    
    opts1 = {}
    opts1.update(kw)
    opts1.setdefault("createParents",True)
    if opts1.pop("createParents"):
        try:
            dirName = opts1["dir"]
        except KeyError:
            pass
        else:
            makedir(dirName)
    l2 = []
    opts1.setdefault("mode","w")
    for k in ("mode","bufsize"):
        if opts1.has_key(k):
            l2.append(opts1[k])
            del opts1[k]
    if opts1.pop("withTime",False):
        opts1["prefix"] = opts1.get("prefix","") + currTimeAsFileName()+'.'
    (fd,name) = mkstemp(*l,**opts1)
    return (os.fdopen(fd,*l2),name)

def makeWorkFile(pathBase,returnPathOnly=True):
    dirName,baseName = os.path.split(pathBase)
    fobj,path = makeTmpFile(suffix='.tmp', prefix=baseName, dir=dirName)
    if returnPathOnly:
        fobj.close()
        return path
    else:
        return (fobj,path)

def strToFile(s,fileName,mode="w",dryRun=False):
    s = str(s)
    if not dryRun:
        out = openCompressed(fileName,mode)
        out.write(s)
        out.close()
    else:
        print "Write to file %s mode %s:" % (fileName,mode)
        print s

def fileToStr(fileName):
    inp = openCompressed(fileName,"r")
    s = inp.read()
    inp.close()
    return s

def fileSync(f):
    """Flush both user and kernel buffers of an open file object.
    @param f file object as returned by open()"""
    f.flush()
    os.fsync(f.fileno())

def stripSfx(s,sep='.'):
    """Remove right-most instance of separator string and everything past it.
    Primary use is to remove the file name 'extension'.
    @return input string s without the suffix or original input if suffix is not found"""
    return s.rsplit(sep,1)[0]

def stripPathSfx(s,sep='.'):
    """Same as stripSfx but looks only inside the os.path.basename(s)"""
    (head,tail) = os.path.split(s)
    return os.path.join(head,stripSfx(tail,sep=sep))

def strFilter(s,allowed):
    allowed = set(allowed)
    return ''.join(( x for x in s if x in allowed ))

_osNameFilter_allowedDef = set(string.ascii_letters + string.digits + '-_+#$~^.,')

def osNameFilter(s,allowed=_osNameFilter_allowedDef,remove=''):
    return strFilter(s,allowed=allowed-set(remove))

def strToFileName(s,remove=''):
    """Rough conversion of artbitrary string with spaces into a file name w/o spaces.
    Possibly a many-to-one operation."""
    return osNameFilter(s.replace(' ','_'),remove=remove)

class SymbolRunsCompressor:
    """Compress all consequitive runs of identical symbol(s) sym that are longer than minLen into minLen."""
    def __init__(self,sym,minLen):
        #assert len(sym) == 1
        self.rex = re.compile('[%s]{%i,}'%(sym,minLen+1))
        self.rep = sym*minLen

    def __call__(self,s):
        """Compress the sequence.
        @param s string or numpy character array
        """
        if isinstance(s,str):
            return re.sub(self.rex,self.rep,s)
        else: #numpy array
            return n.fromstring(re.sub(self.rex,self.rep,s.tostring()),dtype="S1")

class ArrStr:
    """A collection of string methods on numpy string_ arrays to insulate the code from changes in numpy interfaces.
    @todo For newer numpy versions, it should just redirect to numpy.char free functions"""
    @staticmethod
    def upper(s):
        return n.fromstring(s.tostring().upper(),dtype='S1')

def isSamePath(path1,path2):
    paths = [ os.path.abspath(os.path.realpath(p)) for p in (path1,path2) ]
    return paths[0] == paths[1]

def editSymlink(source,link_name):
    """Create the new symlink or change existing symlink.
    @param source What symlink will point to
    @param link_name Path to symlink itself
    """
    if os.path.islink(link_name):
        os.remove(link_name)
    os.symlink(source,link_name)


def openCompressed(filename,mode,compressFormat=None,**kw):
    """Open a filename which can be either compressed or plain file.
    @param compressFormat if None, an attempt will be made to autodetect format 
    (currently by extension, only '.gz' and '.bz2' are recognized); if "none" - open as plain
    file, if "gzip" - open as gzip file."""
    cform = compressFormat
    if cform is None:
        cform = "none"
        if filename.endswith('.gz'):
            cform = "gzip"
        elif filename.endswith('.bz2'):
            cform = "bz2"
    #print "DEBUG: openCompressed(%s,%s,%s)" % (filename,mode,cform)
    k = kw.copy()
    ret = None
    if cform == "gzip":
        if "buffering" in k:
            k["bufsize"] = k["buffering"]
            del k["buffering"]
        ret = openGzip(filename,mode,**kw)
    elif cform == "bz2":
        k.setdefault("buffering",2**20)
        ret = bz2.BZ2File(filename,mode,**kw)
    elif cform == "none":
        k.setdefault("buffering",2**20)
        ret = open(filename,mode,**kw)
    else:
        raise ValueError(compressFormat)
    return ret

def compressFile(inputFile,outputFile=None,compressFormat="gz"):
    """Compress input file in place, like gzip does by default.
    Optionally, move the result into outputFile"""
    if compressFormat == "gz":
        run(["gzip","--force",inputFile])
        if outputFile:
            shutil.move(inputFile+".gz",outputFile)
    else:
        raise ValueError("Unknown compressFormat: %s" % (compressFormat,))


class PopenStdinProxy(object):
    """Proxy class to use in place of Popen.stdin - it makes sure than after stdin.close(), .wait() is called on the subprocess.
    When Python exits, it apparently terminates still running subprocesses (maybe SGE does it - I only observed the issue in
    batch execution. In particular, when we call .close() on a pipe to gzip subprocess writer, it takes some time for gzip
    to finish writing. It our application exits immediately, gzip is killed leaving unfinished output file.
    This proxy class holds a reference to the Popen instance and delegates all attribute lookups to Popen's stdin file object,
    except for the .close() method, in which it first closes the stream, and then waits for the child process to finish.
    This is not ideal solution because Popen.stdin is a real built-in file object, and this proxy class is not. Thus, some
    code that e.g. uses Python C Api can choke on it. Another alternative would be to keep a global list of all Popen
    objects and wait on them in atexit() handler.
    """

    def __init__(self,p):
        self.p = p

    def __getattribute__(self,x):
        if x == "close":
            return object.__getattribute__(self,"close")
        else:
            return getattr(object.__getattribute__(self,"p").stdin,x)

    def close(self):
        p = object.__getattribute__(self,"p")
        p.stdin.close()
        p.wait()
        if p.returncode:
            raise CalledProcessError("",p.returncode)

    def __del__(self):
        self.close()

def openGzip(filename,mode,compresslevel=6):
    compresslevel = int(compresslevel)

    if mode in ("w","wb","a","ab"):
        if mode in ("w","wb"):
            redir = ">"
        elif mode in ("a","ab"):
            redir = ">>"
        #p = Popen("gzip -%s %s %s" % (compresslevel,redir,filename), shell=True, env=os.environ, bufsize=2**16, stdin=PIPE, close_fds=True)
        stdout = open(filename,mode)
        p = Popen(["gzip","-%s" % compresslevel,"-c"], shell=False, env=os.environ, bufsize=2**16, stdin=PIPE, stdout=stdout, close_fds=True)
        stdout.close()
        return PopenStdinProxy(p)

    elif mode in ("r","rb"):
        cmd = ["gzip","-cd",filename]
        if not os.path.isfile(filename):
            raise OSError("Input file does not exist: %s", filename)
        p = Popen(cmd,env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True)
        if p.returncode:
            raise CalledProcessError(str(cmd),p.returncode)
        return p.stdout
    
    else:
        raise ValueError("'openGzip()' - Unsupported 'mode': " + mode)

_ioFilterCodeTpl = """import sys
%s
_iofilter_code=%s
for line in sys.stdin:
    sys.stdout.write(_iofilter_code(line))
"""

def ioFilter(inp,code,mode="stream",lineInitCode=""):
    """Pipe input stream through Python code and return in the output stream.
    @param inp Input file object
    @param code Python code (possibly as multi-line) string that will be 
    executed under a separate Python instance; it should read from standard input 
    and write to standard output
    @param mode if "line", then code accepts one line and returns one or more lines, 
    e.g. "lambda x: x.replace(',','\n')";
    if "stream", then code reads from standard input and writes to standard output.
    @return Input file object connected to the standard output of the filter
    @param lineInitCode if mode == "line", code will be passed to eval() to 
    prepare a callable object, and thus must conform to eval()'s restrictions 
    (the main in this context would be inability to use 'import'). To address 
    this restriction, and optional argument lineInitCode can be inserted before 
    the main loop and can contain module imports or method definitions.
    """
    assert mode in ("stream","line"),"Unknown stream parameter value: %s" % (mode,)
    codeS,codeFileName = makeTmpFile(prefix="iofilt",withTime=True)
    if mode == "stream":
        codeS.write(code)
    elif mode == "line":
        codeS.write(_ioFilterCodeTpl % (lineInitCode,code))
    codeS.close()
    out = Popen([sys.executable,codeFileName],env=os.environ, bufsize=2**16, stdin=inp, stdout=PIPE, close_fds=True).stdout
    #os.remove(codeFileName)
    #print codeFileName
    return out

def strAttributes(o,exclude=tuple(),delim=' | '):
    """Return a string with all attributes names and value for an object 'o'."""
    return delim.join([ name + ' : ' + str(val) for (name,val) in o.__dict__.items()
        if not name in exclude])

def delAttr(o,attr,silent=True):
    try:
        delattr(o,attr)
    except AttributeError:
        if not silent:
            raise


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

def masksToInd(maskVal,maskInd):
    """Convert value[0:someN] and index[0:someN] into value[index]. 
    Assumes N->1 relation between index and value. In particular, that means
    value[i] == 0 && value[j] != 0 && index[i] == index[j] is not allowed.
    No checking is made for the above prerequsite - the caller is responsible."""
    val = numpy.zeros(numpy.max(maskInd)+1,dtype=maskVal.dtype)
    val[maskInd] = maskVal
    return val
    
def whereItems(arr,condition):
    wh = numpy.where(condition)
    return numpy.rec.fromarrays((wh[0],arr[wh]),names='ind,val')

def fromWhereItems(whItems,defVal=0):
    wh =  whItems['ind']
    a = numpy.empty(numpy.max(wh) + 1, dtype = whItems['val'].dtype)
    a[:] = defVal
    a[wh] = whItems['val']
    return a

def logicalAnd(*arrays):
    """Does the same as numpy.logical_and(), but works on a list of arrays of arbitrary length"""
    res = n.logical_and(arrays[0],arrays[1])
    if len(arrays) == 2:
        return res
    for a in arrays[2:]:
        res = n.logical_and(res,a)
    return res

def allSame(seq,key=None,comp=None):
    """Return True if all elements of iterable are equal.
    @param seq Any iterable
    @param key Unary operator that gives the key to compare as key(x) for x in seq. None means identity.
    @param comp Comparison equality binary operator (transitivity holds). None means ==.
    @note Will raise ValueError on empty iterable
    @note It assumes transitivity of equality comparison, and only checks first element with all others
    @note Use case: allSame(arrays,key=len) will check that all arrays have the same length
    @note Use case: allSame(arrays,comp=lambda x,y: numpy.all(x==y))
    """
    if key is None:
        key = lambda x: x
    if comp is None:
        comp = lambda x,y: x == y
    itr = ( key(x) for x in seq )
    try:
        first = itr.next()
    except StopIteration:
        raise ValueError("Need non-empty iterable")
    for x in itr:
        if not comp(first,x):
            return False
    return True
    
def diffFiles(fromFile,toFile,showCommon=False,asText=False):
    """Compare two files line-by-line using difflib module.
    @param fromFile Diff is computed from this file
    @param toFile Diff is computed toward this file
    @param showCommon If False [default], only show lines that changed.
    difflib by default shows identical lines as well.
    @param asText If True, return result as a single multiline string; otherwise
    return as an iterator for lines of output (default).
    @return A generator for diff strings. See docs on difflib for the format
    of output lines.
    @note Currently we use difflib.ndiff()
    """
    import difflib
    def _file_to_lines(f):
        inp = openCompressed(f,"r")
        lines = inp.readlines()
        inp.close()
        return lines
    fromLines = _file_to_lines(fromFile)
    toLines = _file_to_lines(toFile)
    if not showCommon:
        filt = lambda line: not line.startswith('  ')
    else:
        filt = lambda line: True
    outIter = ( line for line in difflib.ndiff(fromLines,toLines) if filt(line) )
    if asText:
        return ''.join(outIter)
    else:
        return outIter

def sameArrays(arrays):
    return allSame(arrays,comp=lambda x,y: n.all(x==y))

def binCount(seq,format="dict"):
    cnt = defdict(int)
    #numpy recarray records will be compared by address (not what we need),
    #so we convert everything numpy to python builtins
    for x in seq:
        cnt[x.item() if isinstance(seq[0],n.generic) else x] += 1
    if format == "dict":
        return cnt
    elif format == "list":
        return sorted(cnt.items())
    else:
        raise ValueError("unknown format value: %s" % format)

def recFromArrays(arrays,names):
    """Does the same as numpy.rec.fromarrays() but does not fail when arrays are themselves recarrays.
    Also, if any of the arrays is a numpy.rec.record scalar, it will be broadcast to the length
    of other arrays (same as in array assignment semantics). Always creates a fresh copy of the data.
    @param arrays sequence of numpy arrays (possibly record arrays themselves), of the same shape
    @param names sequence or a comma delimited string with new field names, one per input array
    @return new record array with concatenated records"""
    if isinstance(names,str):
        names = names.split(",")
    assert len(names) == len(arrays) and len(arrays)>0, "'arrays' and 'names' must have non-zero and equal length"
    arrMax = arrays[n.argmax([sum(arr.shape) for arr in arrays])]
    for arr in arrays:
        assert arr.shape == arrMax.shape or len(arr.shape) <= 1,\
            "Input arrays must have the same shape or be numpy scalars or size 1 arrays"
    dt = [(name,arr.dtype) for (name,arr) in it.izip(names,arrays)]
    out = n.empty(arrMax.shape,dtype=dt)
    for (name,arr) in it.izip(names,arrays):
        out[name] = arr
    return out

def flattenDtype(dtype,withPrefix=False):
    """Given a nested record dtype, return a flat dtype record dtype.
    @param dtype recarray dtype
    @param withPrefix if True, add each compound field name as prefix to the contained field name.
    Example: if dtype is 
    [('param', [('C', '<f8'), ('thresh', '<f4')]), ('val', [('senMean', '<f4'), ('speMean', '<f4')])]
    and withPrefix is False, the return value is
    [('C', '<f8'), ('thresh', '<f4'), ('senMean', '<f4'), ('speMean', '<f4')]
    if withPrefix is True, the return value is
    [('param_C', '<f8'), ('param_thresh', '<f4'), ('val_senMean', '<f4'), ('val_speMean', '<f4')]
    if withPrefix is False and field names are not unique, the result is undefined."""
    flatDescr = []
    def _flattenDescr(pref,descr,flatDescr):
        assert isinstance(descr,list),"'union' dtypes are not supported"
        for el in descr:
            if len(el) == 2: #(name,dtype)
                name = el[0]
                dt = el[1]
            elif len(el) == 3: #(title,name,dtype)
                name = el[1]
                dt = el[2]
            if isinstance(dt,list):
                _flattenDescr(pref+name+'_',dt,flatDescr)
            else:
                assert isinstance(dt,str)
                if withPrefix:
                    pr = pref
                else:
                    pr = ''
                flatDescr.append((pr+name,dt))
    _flattenDescr('',dtype.descr,flatDescr)
    return n.dtype(flatDescr)



def selFieldsArray(arr,fields):
    """Select a list of fields from a record array as a new record array.
    @todo current method creates a full copy; for ajacent fields, that
    can be done as a view by splitting all fields into two sub-dtypes."""
    return recFromArrays([arr[f] for f in fields],names=','.join(fields))

def joinRecArrays(arrays):
    """Same as SQL table join with row number used as a join index.
    @param arrays one or more recarrays, with all fields names unique across all arrays
    @return record array with fields being a union of fields in input arrays
    """
    #If non-unique names are present, the numpy error message will be cryptic,
    #so we do error checking here.
    #dtype.descr might have title before name in dtype.descr, so we need to use
    #accessor attribute to get to names reliably
    if not isUniqueArray(reduce(operator.add,[arr.dtype.names for arr in arrays])):
        raise ValueError,"All names in input record arrays must be unique"
    dt = reduce(operator.add,[arr.dtype.descr for arr in arrays])
    grNames = [ "f%i" % i for i in range(len(arrays)) ]
    #field name idexing is only allowed with a single field argument
    #we use this trick to access a group of ajacent fields as a whole:
    #create the new array with compound fields - each compound field has
    #all fields from one input array, and then return a view with a flattened
    #field list.
    #@todo write a separate method that returns a compound view for selected groups
    #of fields - can be usefull in any application that needs to access (e.g. fast assign)
    #a group of fields at once.
    return recFromArrays(arrays,names=grNames).view(dt)

def fieldDtypeRecArray(arr,name):
    """Return dtype of a field with name 'name' from a record array 'arr'.
    Note: you should not try to determine field's dtype as a single element's
    dtype, such as arr[0][name].dtype, because this will return the dtype
    of a numpy scalar, which for string fields will the length of a given
    value that might be smaller than the field length in a record array."""
    return arr.dtype.fields[name][0]

def permuteObjArray(arr):
    #numpy permutation or shuffle do not work on arrays with 'O' datatypes
    return arr[nrnd.permutation(n.arange(len(arr),dtype=int))]

def groupPairs(data,keyField=0):
    """Create a dict(key->list of vals) out of a sequence of pairs (key,val).
    @param data sequence of pairs
    @param keyField index of field to use as key (the remaining field is used as val)"""
    x = defdict(list)
    assert keyField in (0,1)
    valField = 1 - keyField
    for rec in data:
        key = rec[keyField]
        val = rec[valField]
        x[key].append(val)
    # important to convert to regular dictionary otherwise we get
    # silent insertion of default elements on access to non-existing keys
    return dict(x)

def groupTuples(data,keyField=0):
    """Create a dict(key->list of tuples) out of a sequence of tuples.
    Similar to creating a non-unique database index or sorting by a given field, 
    with the exception that the data is copied inside the index.
    @param data sequence of tuples
    @param keyField index of field to use as key (all fields are used as val)"""
    x = defdict(list)
    for rec in data:
        x[rec[keyField]].append(rec)
    # important to convert to regular dictionary otherwise we get
    # silent insertion of default elements on access to non-existing keys
    return dict(x)

def indexTuples(data,keyField=0):
    """Create a dict(key->tuple) out of a sequence of tuples.
    Similar to creating a unique database index. 
    A ValueError exception will be raised if records are not unique.
    @param data sequence of tuples
    @param keyField index of field to use as key (all fields are used as val)"""
    x = dict()
    #in order to test for uniqueness we need the length of the input,
    #but len() will not work if the input is an iterator, so we count
    for (i_rec,rec) in enumerate(data):
        x[rec[keyField]] = rec
    if len(x) <= i_rec:
        raise ValueError("Non-unique items in input data")
    return x

def recFromRecords(recs,dtype=None):
    """Does the same as numpy.rec.fromrecords() but does not fail when records are recarrays with 'O' fields.
    This works around a bug in numpy 1.3.
    Always creates a fresh copy of the data.
    @param recs sequence of numpy rec array records or tuples.
    @param dtype dtype for the new array, if None - will be taken from recs[0], 
    which in that case must be non-empty list of numpy recs
    @return new record array with concatenated records
    @todo probably use in groupRecArray()"""
    n_recs = len(recs)
    if dtype is None and n_recs == 0:
        raise ValueError("Need dtype specified when an empty sequence of records is given")
    out = n.empty(n_recs,dtype=recs[0].dtype if dtype is None else dtype)
    for (i_rec,rec) in enumerate(recs):
        out[i_rec] = tuple(rec)
    return out

def has_numpy_bug_object_cast_rec_arr(arr):
    """Check if a bug when numpy cannot convert a rec array with 'O' fields to itself with dtype argument"""
    # we need at least one record to check
    if len(arr) == 0:
        return False
    else:
        # bug happens when we build a list of records and try to convert it to array
        recs = [ arr[0] ]
        try:
            n.asarray(recs,dtype=arr.dtype)
        except ValueError:
            return True
        else:
            return False

def numpyToScalarFunc(x):
    """Return a function object that will be either numpy.asscalar or identity depending on type of x"""
    try:
        n.asscalar(x)
        _to_sc = n.asscalar
    except AttributeError:
        _to_sc = lambda x: x
    return _to_sc

def groupRecArray(arr,keyField):
    """Create a dict(key->record array) out of record array.
    Similar to creating a non-unique database index or sorting by a given field, 
    with the exception that the data is copied inside the index.
    @param arr Numpy record array
    @param keyField name of field to create the key from
    @return dict that for each unique value of keyField contains a numpy record 
    array with all records that have this key value.
    @post Grouping is stable - the original order of records within each group is preserved.
    @bug This will not work if there are fields of type "O" - seems to be a numpy bug"""
    m = defdict(list)
    if len(arr):
        _tokey = numpyToScalarFunc(arr[0][keyField])
        if has_numpy_bug_object_cast_rec_arr(arr):
            for rec in arr:
                #dereferenced recarray field scalars are compared by address in dict (because they are mutable?), 
                #hence asscalar (which calls .item()) to get immutable Python scalar
                m[_tokey(rec[keyField])].append(tuple(rec)) #have numpy bug, convert rec to tuple
        else:
            for rec in arr:
                #dereferenced recarray field scalars are compared by address in dict (because they are mutable?), 
                #hence .item() to get immutable Python scalar
                m[_tokey(rec[keyField])].append(rec)
        for key in m:
            m[key] = n.asarray(m[key],dtype=arr.dtype)
        # important to convert to regular dictionary otherwise we get
        # silent insertion of default elements on access to non-existing keys
    return dict(m)

def indexRecArray(arr,keyField):
    """Create a dict(key->record of array) out of record array.
    Similar to creating a unique database index. 
    A reference to a single array record is stored as a value field in each item of the index.
    A ValueError exception will be raised if records are not unique.
    @param arr Numpy record array
    @param keyField name of field to create the key from
    @return dict that for each unique value of keyField contains a rerefernce to a single numpy record."""
    m = {}
    for rec in arr:
        #dereferenced recarray field scalars are compared by address in dict (because they are mutable?), 
        #hence .item() to get immutable Python scalar
        m[rec[keyField].item()] = rec
    if len(m) < len(arr):
        raise ValueError("Non-unique items in input array")
    return m

def countRecArray(arr,keyFields,format="dict"):
    return binCount(selFieldsArray(arr,fields=keyFields),format=format)

def isUniqueArray(arr):
    return len(n.unique(arr)) == len(n.ravel(arr))


def saveRecArrayAsCsv(arr,out,withHeader=True,sep=','):
    outDt = flattenDtype(arr.dtype,withPrefix=True)
    outArr = arr.view(outDt)
    closeOut = False
    if isinstance(out,str):
        out = open(out,'w')
        closeOut = True
    if withHeader:
        out.write(sep.join(outDt.names)+'\n')
    for rec in outArr:
        out.write(sep.join([ "%s" % (x,) for x in rec ])+'\n')
    if closeOut:
        out.close()

def rndRound(x):
    """Round real non-negative scalar to the nearest integer in a random way.
    The closer x to its ceil(), the larger the probability for
    x to be rounded to ceil().
    This is useful when we need to pick a multiple times count of random samples
    which is a fraction of something, such that the total is on average the same.
    @return rounded value of type int
    """
    assert x >= 0
    b = n.floor(x)
    c = x - b
    p = nrnd.random()
    return b + (c > p)

class SubSamplerUniRandomEnd:
    """Uniform random [0,rnd_length] subsampler where rnd_length is in [minLen,maxLen]"""

    def __init__(self,minLen,maxLen):
        assert minLen > 0 and maxLen >= minLen
        self.minLen = minLen
        self.maxLen = maxLen

    def __call__(self,samp):
        """Return subsequence [0,random).
        We always take from the beginning rather than from a random start,
        because when subsampling short taxa with concatenation, this gives 
        a better chance of not hitting spacer junction."""
        sampLen = len(samp)
        return samp[0:nrnd.random_integers(min(sampLen,self.minLen),min(sampLen,self.maxLen))]

class SubSamplerRandomStart:
    """random [rnd_start,rnd_start+rnd_length] subsampler where rnd_length is in [minLen,maxLen]"""

    def __init__(self,minLen,maxLen=None):
        if maxLen is None:
            maxLen = minLen
        assert minLen > 0 and maxLen >= minLen
        self.minLen = minLen
        self.maxLen = maxLen

    def __call__(self,samp):
        """Return subsequence [random,random)."""
        sampLen = len(samp)
        fragLen = nrnd.random_integers(min(sampLen,self.minLen),min(sampLen,self.maxLen))
        fragStart = nrnd.random_integers(0,sampLen-fragLen)
        return samp[fragStart:fragStart+fragLen]


def getRectStencil(arr,center,halfSize):
    """Return a stencil from Numpy 2D array.
    @todo looks like numpy.ix_ can do that and more"""
    center = n.asarray(center)
    ind = n.clip((center-halfSize,center+halfSize+1),(0,0),n.asarray(arr.shape)-1)
    sel = arr[ind[0][0]:ind[1][0],ind[0][1]:ind[1][1]]
    #print center, halfSize, ind, sel
    return sel,ind

class ArrayAppender:
    """Wrapper for Numpy 1D array that emulates list.append() behaviour"""
    def __init__(self,arr):
        self.arr = arr
        self.num = 0

    def getMemory(self):
        """Return entire internal array"""
        return self.arr

    def getData(self):
        """Return only used portion of the internal array"""
        return self.arr[:self.num]

    def __len__(self):
        """Return length of getData() result"""
        return self.num

    def nextItem(self):
        """Advance to the next unused element, resizing (2x) the internal array when necessary
        @return tuple (internal array, index of the next element)"""
        inext = self.num
        self.num += 1
        arr = self.arr
        if inext >= len(arr):
            arr = n.append(arr,n.zeros(len(arr),arr.dtype))
            self.arr = arr
        return (arr,inext)

    def nextElem(self):
        """Like nextItem(), but return the next element itself.
        That only makes sense when the array has mutable elements (e.g. is a record array)"""
        arr,inext = self.nextItem()
        return arr[inext]


def enumNames(obj):
    """Return a tuple with all attribute names of an obj not starting with '_'
    @param obj instance or class"""
    return tuple([ name for name in dir(obj) if not name.startswith('_') ])

def enumNamesToValues(names,obj):
    """Convert a list of attribute names to attribute values for a given object.
    This is designed to process mnemonic enum constants into values,
    assuming enums are implemented as class or instance attributes.
    Example: if enums are defined as:
    class ENUM_FILE_MODE:
        WRITE = 0x01
        READ = 0x02
        BINARY = 0x10
        TEXT = 0x20
    and out program accepts an option such as "--run-mode WRITE,BINARY",
    the call to enumNamesToValues(opt.runMode,ENUM_FILE_MODE)
    will return (0x01,0x10) while also checking that each name is a valid
    integer attribute and does not start with undescore. ValueError or 
    AttributeError will be raised if these conditions are not met.
    In a typical case when the values are bitmasks that can be ORed,
    use listOr(enumNamesToValues(names,obj))
    @param sequence with names or a comma separated string of names
    @param obj object from which to extract the attribute values - can be instance or class
    @return tuple of values"""
    oneName = False
    if isinstance(names,str):
        names = names.split(',')
    names = [ name.strip() for name in names ]
    vals = []
    for name in names:
        if name.startswith('_'):
            raise ValueError("Name cannot start with underscore: %s" % (name,))
        val = getattr(obj,name)
        if not isinstance(val,int):
            raise ValueError("Enum value must have integer type: %s = %s" % (name,val))
        vals.append(val)
    return tuple(vals)

def listOr(l):
    return reduce(operator.or_,l)

def binIntToStr(x,digits=8):
    """Return a string representation of x as a binary encoded number.
    Taken from some discussion list."""
    return ''.join(x & (1 << i) and '1' or '0' for i in range(digits-1,-1,-1))

def floatsToPcntStr(arr):
    return [ "%.2f" % (x*100) for x in arr ]

def floatToPcntStr(x):
    return "%.2f" % (x*100)

def dictToMaskedArray(d,maxInd=None,dtype='O'):
    """Create masked array out of dict
    @param d dict(index->value) where index is integer type
    @param maxInd maxiumum possible index value plus one, if None, obtain from dict keys
    @param dtype dtype of the result
    @return Numpy masked array with values not present in the original dict masked"""
    if maxInd is None:
        maxInd = int(max(d.keys()))+1
    ind = n.zeros(maxInd,dtype=dtype)
    mask = n.zeros(maxInd,dtype=bool)
    mask[:] = True
    for (i,v) in d.items():
        i = int(i)
        ind[i] = v
        mask[i] = False
    return nma.masked_array(ind,mask=mask)

def izipCount(l):
    """Return iterator of pairs (0,l[0]),(1,l[1]),... where l is an input iterable"""
    return it.izip(it.count(),l)

def runsOfOnesArray(bits):
    """Return Nx2 array where each row specifies the range for a run of 1 (or True) in 'bits' input array"""
    # This is taken from
    # http://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array
    # This might be made slightly more memory efficient if the alternative recipe is used from the bottom of 
    # the referenced page:
    # >>> from numpy import array, arange
    # >>> b = array([0,0,0,1,1,1,0,0,0,1,1,1,1,0,0], dtype=bool)
    # >>> sw = (b[:-1] ^ b[1:]); print sw
    # [False False  True False False  True False False  True False False False
    #           True False]
    # >>> isw = arange(len(sw))[sw]; print isw
    # [ 2  5  8 12]
    # >>> lens = isw[1::2] - isw[::2]; print lens
    # [3 4]
    # But that will require some processing depending on the start element.

    # make sure all runs of ones are well-bounded
    bounded = numpy.hstack(([0], bits, [0]))
    # get 1 at run starts and -1 at run ends
    difs = numpy.diff(bounded)
    run_starts, = numpy.where(difs > 0)
    run_ends, = numpy.where(difs < 0)
    # go back to offsets before padding
    #run_starts -= 1
    #run_ends -= 1
    return n.column_stack([run_starts,run_ends])

def printAliSeqs(seqs,lineLen,out,seqNames=None,emptySymb=' '):
    """Print two aligned strings side-by-side"""
    maxLen = max((len(s) for s in seqs))
    if seqNames is None:
        seqNames = ['']*len(seqs)
    seqs = [ s.ljust(maxLen,emptySymb) for s in seqs ]
    maxNameLen = max((len(s) for s in seqNames))
    for x in range(0,maxLen,lineLen):
        for (s,sName) in it.izip(seqs,seqNames):
            out.write(sName.ljust(maxNameLen+1))
            out.write(s[x:x+lineLen])
            out.write("\n")
        out.write("\n")

def openCsv(csvFile,mode,factory=None,*l,**kw):
    import csv
    if factory is None:
        if mode.startswith("r"):
            factory = csv.reader
        else:
            factory = csv.writer
    csvFileStream = None
    if isinstance(csvFile,str):
        csvClose = True
        csvFileStream = openCompressed(csvFile,mode)
        csvFile = factory(csvFileStream,*l,**kw)
    else:
        csvClose = False
        #Here we need to figure out if csvFile is just a file stream, or
        #a CSV reader/writer already. Unfortunately, csv module does not specify
        #any common base class for CSV objects, so we have to rely on
        #tests for attribute presence
        if not (hasattr(csvFile,"dialect") \
                and (hasattr(csvFile,"line_num") or hasattr(csvFile,"writerow"))):
            #this is NOT a result of calling csv.reader() or compatible interface,
            #so we assume it to be a file stream object, and call factory()
            #on it to create a CSV object
            csvFile = factory(csvFile,*l,**kw)
    return Struct(csvClose=csvClose,csvFile=csvFile,csvFileStream=csvFileStream)

class ObjCache(object):
    """Class that caches immutable objects based on identity.
    Use case: we are reading N:1 relation from a two-column csv file
    where both columns are represented by fairly long string IDs.
    E.g. read id -> contig id.
    If we use these string as read, each one will allocate a separate 
    object. This is wastefull if there is a lot of repeated values
    (e.g. for contig id if there are many reads in the same contig).
    This object keeps a hash of objects and always returns a reference to
    an existing one."""

    def __init__(self):
        self.store = dict()

    def __call__(self,o):
        return self.store.setdefault(o,o)

    def free(self,o):
        del self.store[o]

    def clear(self):
        self.store.clear()

## Global instance of object cache
objCache = ObjCache()

def join_sorted_records_right_check(iter1,iter2,key1,key2):
    """Join two sorted iterator ranges checking that that there is always a match on the right.
    @pre input records are sorted, each left record has match on the right.
    @param iter1 iterator of record tuples sorted by key1
    @param iter2 iterator of record tuples sorted by key2
    @param key1 Identifies key in iter1 records. Int means index, otherwise must be functor key(record)
    @param key2 Identifies key in iter2 records. Int means index, otherwise must be functor key(record)
    @return iterator of (record1,record2), where key(record1)==key(record2). Other records are skipped.
    """
    #note: more general kinds of merges can be implemented by using itertools.group(heapq.merge()), 
    #but this specific case is simple enough.
    if isinstance(key1,int):
        key1 = lambda rec,_key=key1: rec[_key]
    if isinstance(key2,int):
        key2 = lambda rec,_key=key2: rec[_key]
    for rec1 in iter1:
        for rec2 in iter2:
            if key1(rec1) == key2(rec2):
                yield (rec1,rec2)
                break
        #this is a rarely used 'else' at the end of the loop keyword;
        #it will be reached if 'break' is never executed
        else:
            raise KeyError(key1(rec1))

def basedir(path):
    """Return right-most component of the directory path.
    This gives identical results both for /some/base and /some/base/"""
    dirn,basen = os.path.split(path)
    if not basen:
        basen = os.path.basename(dirn)
    return basen

def split_path(path):
    """Split path in a list of components similar to str.split().
    Unlike os.path.split(), this will split 'a/b/c' into ['a','b','c'].
    Otherwise, semantics is the same on various corner cases if
    the last returned component as basename and all the leading ones -
    as dirnames. In particular, at least two components are always returned,
    with one of them possibly empty.

    os.path.join(split_path(path)) == path"""

    parts = []
    pathn = path
    while True:
        dirn,basen = os.path.split(pathn)
        if dirn == pathn:
            if len(parts) == 0:
                parts.insert(0,basen)
            parts.insert(0,dirn)
            break
        parts.insert(0,basen)
        if dirn:
            pathn = dirn
        else:
            if len(parts) == 1:
                parts.insert(0,dirn)
            break
    return tuple(parts)

def urljoin_path(base,url):
    import urlparse
    #urlparse.urljoin is weird: 
    #In [6]: urljoin(urljoin("/","static/vicvb"),"jbrowse")
    #Out[6]: '/static/jbrowse'
    if not base.endswith(":"):
        if not base.endswith("/"):
            base += "/"
    return urlparse.urljoin(base,url)

def to_url_params(params):
    """You might have to pass OrderedDict if the order of parameters
    is important. Alternatively, urllib.quote_plus can be applied
    directly to a string."""
    import urllib
    return urllib.urlencode(params)

def add_to_path(dir,var="PATH",prepend=False,env=None):
    """Add a directory to the PATH environment variable"""
    dir = str(dir)
    if env is None:
        env = os.environ
    if var in env:
        if prepend:
            first = dir
            second = env[var]
        else:
            first = env[var]
            second = dir
        env[var] = os.pathsep.join((first,second))
    else:
        env[var] = dir

def tar_check_safety(tar):
    
    def _tar_info_str(tarinfo):
        return " ; ".join([ "%s : %s" % item for \
                item in sorted(tarinfo.__dict__.items()) \
                if not item[0].startswith('_') \
                and not item[0] == "buf" ])
    
    def _err_msg(tarinfo,msg):
        return "Archive failed safety check - "+\
                    msg+": %s" % (_tar_info_str(tarinfo),)

    for tarinfo in tar:
        if os.path.isabs(tarinfo.name) or \
                os.path.isabs(os.path.normpath(tarinfo.name)):
            raise ValueError(_err_msg(tarinfo,
            "Absolute file name detected"))
        elif ".." in tarinfo.name or ".." \
                in os.path.normpath(tarinfo.name):
            raise ValueError(_err_msg(tarinfo,
                    "Upper directory reference is detected"))
        elif not (tarinfo.isreg() or tarinfo.isdir()):
            #e.g. if archive was artificially manipulated to contain 
            #first A/B where B is a symlink to ../../something,
            #and then A/B/C, then C might be created as ../../something/C 
            #(my guess).
            raise ValueError(_err_msg(tarinfo,
                    "Non-regular files or dirs can lead to exploits"))

def tar_extractall_safe(archive,path=None):
    if path is None:
        path = os.getcwd()
    tar = tarfile.open(archive, "r") #will auto-detect compression
    try:
        tar_check_safety(tar)
        tar.extractall(path=path)
    finally:
        tar.close()


def tar_extractall_safe_single_dir(archive,path=None):
    tar_extractall_safe(archive=archive,path=path)
    subdirs = list(os.listdir(path))
    assert len(subdirs) == 1,\
            "Expected a single directory in archive %s" \
            % (path,)
    return os.path.join(path,subdirs[0])

def load_config_json(config_file,default={}):
    if os.path.exists(config_file):
        with open(config_file,'r') as f:
            return json.load(f)
    else:
        return default

def save_config_json(config,config_file):
    with open(config_file,'w') as f:
        json.dump(config,f)

def none_from_str(s):
    if s is not None:
        if s == "None":
            return None
    return s

def get_current_script():
    return os.path.abspath(sys.argv[0])

def get_current_python():
    return sys.executable

def copytree_ext(src, dst, symlinks=False, ignore=None, copy_stats=True):
    """Recursively copy a directory tree.
    This is a modified clone of shutil.copytree() from Python 2.6.
    Extra arguments allow controlling whether time should be copied.
    For example, when copying to a scratch file system where files
    are automatically deleted based on age, you do not want to have
    the modification time copied.

    The destination directory must not already exist.
    If exception(s) occur, an Error is raised with a list of reasons.

    If the optional symlinks flag is true, symbolic links in the
    source tree result in symbolic links in the destination tree; if
    it is false, the contents of the files pointed to by symbolic
    links are copied.

    The optional ignore argument is a callable. If given, it
    is called with the `src` parameter, which is the directory
    being visited by copytree(), and `names` which is the list of
    `src` contents, as returned by os.listdir():

        callable(src, names) -> ignored_names

    Since copytree() is called recursively, the callable will be
    called once for each directory that is copied. It returns a
    list of names relative to the `src` directory that should
    not be copied.

    The optional copy_stats selects whether to copy time.

    XXX Consider this example code rather than the ultimate tool.

    """
    names = os.listdir(src)
    if ignore is not None:
        ignored_names = ignore(src, names)
    else:
        ignored_names = set()

    os.makedirs(dst)
    errors = []
    for name in names:
        if name in ignored_names:
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree_ext(srcname, dstname, symlinks, ignore)
            else:
                if copy_stats:
                    shutil.copy2(srcname, dstname)
                else:
                    shutil.copy(srcname,dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error), why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except Error, err:
            errors.extend(err.args[0])
    if copy_stats:
        try:
            shutil.copystat(src, dst)
        except OSError, why:
            if WindowsError is not None and isinstance(why, WindowsError):
                # Copying file access times may fail on Windows
                pass
            else:
                errors.extend((src, dst, str(why)))
    if errors:
        raise Error, errors

class null_file:
    """File-like object that does nothing"""
    
    def write(self,s):
        pass

    def close(self):
        pass

#will have to become isinstance(x,str) in Python 3
def is_string(x):
    return isinstance(x,basestring)

def cmd_self(wrapper=""):
    """Return command line and components of command line that 
    would run current Python script.
    @param wrapper [""] is a wrapper script that should set up
    necessary environment variables and then exec its arguments"""
    exe = get_current_python()
    script = get_current_script()
    return (
            {
            "cmd":"{} {} {}".format(wrapper,
                exe,
                script),
            "wrapper":wrapper,
            "exe":exe,
            "script":script,
            "files":[wrapper,exe,script] if wrapper else [exe,script]
            }
           )

def cmd_python(wrapper=""):
    """Return command line and components of command line that 
    would run current Python.
    @param wrapper [""] is a wrapper script that should set up
    necessary environment variables and then exec its arguments"""
    exe = get_current_python()
    return (
            {
            "cmd":"{} {}".format(wrapper,
                exe),
            "wrapper":wrapper,
            "exe":exe,
            "files":[wrapper,exe] if wrapper else [exe]
            }
           )


class PathHasher:
    """Converts file paths into a deep directory hierarchy.
    Use this when you need to create so many files 
    that placing them in a single directory will slow
    down file system meta operations such as search.
    """

    param_file = ".path_hasher.json"

    @classmethod
    def is_instance_dir(klass,root):
        return os.path.isfile(os.path.join(root,klass.param_file))

    def __init__(self,root,n_subdirs=2**10,mode="w"):
        """Ctor.
        @param root Top level directory of the hierarchy.
        When converting input paths, it will be assumed that
        they are given relative to this root.
        @param n_subdirs maximum number of internmediate 
        subdirectories per parent directory
        @param mode {"w","r"}["w"] default access mode. If "w",
        generate new hierarchy; If "r" - expect an existing 
        hierarchy and read the parameters saved under "root".
        """
        import math
        #depth is currently fixed at 1
        #the depth of existing hierachy could be found
        #out from param file by future versions of the code
        param_file = os.path.join(root,self.param_file)
        if mode != "w":
            assert os.path.isfile(param_file)
        if os.path.isfile(param_file):
            opt = load_config_json(param_file)
            self.opt = opt
        else:
            opt = {}
            opt["depth"] = 1
            opt["n_subdirs"] = n_subdirs
            self.opt = opt
            makedir(root)
            save_config_json(opt,param_file)
        self.mode = mode
        n_digits = int(math.ceil(math.log10(opt["n_subdirs"])))
        self.dir_hash_tpl = "ph_{:0%i}" % (n_digits,)
        self.root = root

    def path(self,name,mode=None):
        """Convert a file name into the actual path.
        @param name file name relative to the root directory
        of this object; name can be a multicomponent relative
        path
        @param mode {"w","r"}[self.mode] If "w", intermediate
        directory levels will be created; If "r", intermediate
        directory levels must exist and assertion will be made
        @return join(root,intermediate levels,name)
        
        Examples: 
        a -> root/1/a
        a/b/c -> root/1/a/b/c
        """
        assert not os.path.isabs(name),\
                "Path must be relative (to root of this object): {}".format(name)
        if mode is None:
            mode = self.mode
        #name like "a/b/c" will be hashed only
        #by first component
        first_comp = split_path(name)[0]
        if not first_comp:
            first_comp = os.path.basename(name)
        assert first_comp, "Empty file name is not allowed"
        full_hash = self._gen_bucket_int(hash(first_comp),mode=mode)
        return os.path.join(full_hash,name)

    def buckets(self):
        #when implementing for larger depths, use os.walk and control the depth value
        for sub in os.listdir(self.root):
            if sub.startswith("ph_"):
                psub = os.path.join(self.root,sub)
                if os.path.isdir(psub):
                    yield psub

    def _gen_bucket_int(self,i,mode):
        dir_hash = self.dir_hash_tpl.format( i % self.opt["n_subdirs"])
        full_hash = os.path.join(self.root,dir_hash)
        if mode == "w":
            if not os.path.exists(full_hash):
                #other process could have created full_hash
                #in the meantime, so as long as full_hash exists
                #we ignore any errors
                try:
                    os.makedirs(full_hash)
                except:
                    if not os.path.isdir(full_hash):
                        raise
        else:
            assert os.path.isdir(full_hash),"Bucket directory does not exist: {}".format(full_hash)
        return full_hash

    def random_bucket(self):
        return self._gen_bucket_int(random.randrange(self.opt["n_subdirs"]),mode="w")
    
    def listdir(self):
        for bucket in self.buckets():
            for name in os.listdir(bucket):
                yield os.path.join(bucket,name)

    def glob(self,pattern):
        for bucket in self.buckets():
            for path in glob.iglob(os.path.join(bucket,pattern)):
                yield path

    def walk(self,*l,**kw):
        for bucket in self.buckets():
            for v in os.walk(bucket,*l,**kw):
                yield v

    def __call__(self,name,mode=None):
        return self.path(name,mode=mode)

    def mkdtemp(self,suffix="",prefix=""):
        """Create a temporary directory.
        Directories are spread evenly among the bins.
        @post If ret is a returned value, then 
        self.path(basename(ret)) != ret.
        In other words, you will not be able to map
        the returned file name back to its bucket using
        self.path() method. You should store and use the full
        returned path. The reason behind this is because the
        bucket is determined by applying a hashing function to
        the file name, but the unique file name can be only reliably
        generated without race conditions once the bucket is
        known. To work around this, we first generate a random
        bucket, and then generate
        a unique temporary directory name in the selected bucket.
        @return directory path"""

        return tempfile.mkdtemp(suffix=suffix,
                    prefix=prefix,
                    dir=self.random_bucket())

