### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from types import *
import re
import string, os, sys
pjoin = os.path.join
from time import sleep
import time
from copy import copy, deepcopy
from cPickle import dump, load, Pickler, Unpickler
from cStringIO import StringIO
import numpy
n = numpy
nma = numpy.ma
import numpy.random as nrnd
from tempfile import mkstemp
from textwrap import dedent
from collections import defaultdict as defdict
import operator
import itertools as it
import datetime
import pdb

from MGT.Options import Struct
from MGT.Config import options, Options

from MGT.Run import *

from MGT.Bits.Unique import unique, uniquePick

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
    out = openCompressed(fileName,'w',**kw)
    dump(obj,out,-1)
    out.close()
    
def loadObj(fileName,**kw):
    inp = openCompressed(fileName,'rb',**kw)
    ret = load(inp)
    inp.close()
    return ret

def PickleReader(inp,**kw):
    if isinstance(inp,str):
        closeInp = True
        inp = openCompressed(inp,'rb',**kw)
    else:
        closeInp = False
    pkl = Unpickler(inp)
    while True:
        try:
            x = pkl.load()
        except EOFError:
            break
        yield x
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
        if self.closeOut:
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


def strToFile(s,fileName,mode="w",dryRun=False):
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
    """Compress all consequitive runs of identical symbol sym that are longer than minLen into minLen."""
    def __init__(self,sym,minLen):
        assert len(sym) == 1
        self.rex = re.compile('%s{%i,}'%(sym,minLen+1))
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

def openCompressed(filename,mode,compressFormat=None,**kw):
    """Open a filename which can be either compressed or plain file.
    @param compressFormat if None, an attempt will be made to autodetect format 
    (currently by extension, only '.gz' is recognized); if "none" - open as plain
    file, if "gzip" - open as gzip file."""
    cform = compressFormat
    if cform is None:
        cform = "none"
        if filename.endswith('.gz'):
            cform = "gzip"
    #print "DEBUG: openCompressed(%s,%s,%s)" % (filename,mode,cform)
    if cform == "gzip":
        return openGzip(filename,mode,**kw)
    elif cform == "none":
        return open(filename,mode,2**20,**kw)
    else:
        raise ValueError(compressFormat)

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
        p = Popen("gzip -%s %s %s" % (compresslevel,redir,filename), shell=True, env=os.environ, bufsize=2**16, stdin=PIPE, close_fds=True)
        return PopenStdinProxy(p)

    elif mode in ("r","rb"):
        return Popen(("gzip -cd %s" % (filename,)).split(),env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
    
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
    @param code Python code (possibly as multi-line) string that will be executed under a separate Python instance,
    it should read from standard input and write to standard output
    @param mode if "line", then code accepts one line and returns one or more lines, e.g. "lambda x: x.replace(',','\n')";
    if "stream", then code reads from standard input and writes to standard output.
    @return Input file object connected to the standard output of the filter
    @param lineInitCode if mode == "line", code will be passed to eval() to prepare a callable object, and thus must conform
    to eval()'s restrictions (the main in this context would be inability to use 'import'). Therefore, and optional 
    argument lineInitCode will be inserted before the main loop and can contain module imports or method definitions.
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


def groupRecArray(arr,keyField):
    """Create a dict(key->record array) out of record array.
    Similar to creating a non-unique database index or sorting by a given field, 
    with the exception that the data is copied inside the index.
    @param arr Numpy record array
    @param keyField name of field to create the key from
    @return dict that for each unique value of keyField contains a numpy record 
    array with all records that have this key value.
    @postcondition Grouping is stable - the original order of records within each group is preserved."""
    m = defdict(list)
    for rec in arr:
        #dereferenced recarray field scalars are compared by address in dict (because they are mutable?), 
        #hence .item() to get immutable Python scalar
        m[rec[keyField].item()].append(rec)
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
    return len(n.unique1d(arr)) == len(n.ravel(arr))


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
        """Return subsequence [random,random).
        We always take from the beginning rather than from a random start,
        because when subsampling short taxa with concatenation, this gives 
        a better chance of not hitting spacer junction."""
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

