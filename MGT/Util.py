from types import *
import re
import string, os
from time import sleep
import time
from copy import copy, deepcopy
from cPickle import dump, load
from cStringIO import StringIO
import numpy
n = numpy
import numpy.random as nrnd
from tempfile import mkstemp
from textwrap import dedent
from collections import defaultdict as defdict
import operator
import itertools as it

from MGT.Options import Struct
from MGT.Config import options, Options

from MGT.Run import *

from MGT.Bits.Unique import unique, uniquePick

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

def makeFilePath(fileName):
    """Assume that the argument is a file name and make all directories that are part of it"""
    dirName = os.path.dirname(fileName)
    if dirName not in ("","."):
        makedir(dirName)

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

def stripSfx(s,sep='.'):
    """Remove right-most instance of separator string and everything past it.
    Primary use is to remove the file name 'extension'.
    @return input string s without the suffix or original input if suffix is not found"""
    return s.rsplit(sep,1)[0]


def strFilter(s,allowed):
    allowed = set(allowed)
    return ''.join(( x for x in s if x in allowed ))

_osNameFilter_allowedDef = set(string.ascii_letters + string.digits)

def osNameFilter(s,allowed=_osNameFilter_allowedDef):
    return strFilter(s,allowed=allowed)

class SymbolRunsCompressor:

    def __init__(self,sym,minLen):
        assert len(sym) == 1
        self.rex = re.compile('%s{%i,}'%(sym,minLen+1))
        self.rep = sym*minLen

    def __call__(self,s):
        return re.sub(self.rex,self.rep,s)

def isSamePath(path1,path2):
    paths = [ os.path.abspath(os.path.realpath(p)) for p in (path1,path2) ]
    return paths[0] == paths[1]

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
        assert arr.shape == arrMax.shape or len(arr.shape) == 0,\
            "Input arrays must have the same shape or be numpy scalars"
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
    return x


def groupRecArray(arr,keyField):
    """Create a dict(key->record array) out of record array.
    Similar to creating a non-unique database index or sorting by a given field, 
    with the exception that the data is copied inside the index.
    @param arr Numpy record array
    @param keyField name of field to create the key from
    @return dict that for each unique value of keyField contains a numpy record 
    array with all records that have this key value."""
    m = defdict(list)
    for rec in arr:
        #recarray records are compared by address in dict, for some reason, hence .item()
        m[rec[keyField].item()].append(rec)
    for key in m:
        m[key] = n.asarray(m[key],dtype=arr.dtype)
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

