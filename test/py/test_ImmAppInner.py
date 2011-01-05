"""Test internal behaviour of Glimmer ICM implementation.
We see if the average absolute score values depend on the
length of training sequence.
"""
from MGT.ImmApp import *
from MGT.SeqDbFasta import *

def randomSeq(length):
    # nseq samples, with 1 experiment in each returns nseq x 4 array 
    # with just one '1' in each row, 
    # then argmax along rows gives the winning nucleotide index
    # (seems there is a bug in numpy 1.3.0:
    # when length is int array element (nseq.dtype exists and is 'int64'),
    # nrnd.multinomial(1,p,length) returns a flat array (length,) instead of
    # (length,len(p))
    alph = n.fromstring("ATCG",dtype='S1')
    p = [0.25]*4
    return alph[nrnd.multinomial(1,p,int(length)).argmax(1)]

def genRandomSeqDbIid(path,lengths):
    """Generate random reference sequences as SeqDbFasta.
    Each generated sequence is generated independently"""
    seqDb = SeqDbFasta(path,save=True)
    #assert len(set(lengths)) == len(lengths),"Reference lengths must be unique"
    seqLen = {}
    for (iIcm,length) in enumerate(lengths):
        idIcm = iIcm
        writer = seqDb.fastaWriter(idIcm)
        seq = randomSeq(length)
        writer.record(str(idIcm),seq)
        writer.close()
        seqLen[idIcm] = length
    seqDb.opt.seqLen = seqLen
    seqDb.save()

def genRandomSeqDbCat(path,length,nCopies):
    """Generate random reference sequences as SeqDbFasta.
    The first sequence is generated as random, and each next
    one is a concatenation of N copies of the first."""
    seqDb = SeqDbFasta(path,save=True)
    seqBase = randomSeq(length)
    #assert len(set(lengths)) == len(lengths),"Reference lengths must be unique"
    seqLen = {}
    for iIcm in xrange(nCopies):
        idIcm = iIcm
        writer = seqDb.fastaWriter(idIcm)
        seq = n.tile(seqBase,iIcm+1)
        writer.record(str(idIcm),seq)
        writer.close()
        seqLen[idIcm] = len(seq)
    seqDb.opt.seqLen = seqLen
    seqDb.save()

def genRandomQueryIid(path,lengths):
    writer = FastaWriter(path)
    for iSeq,length in enumerate(lengths):
        seq = randomSeq(length)
        writer.record("%000i_%i" % (iSeq,length),seq)

def genRandomQueryZeroMarkov(path,lengths):
    writer = FastaWriter(path)
    seqBase = randomSeq(max(lengths))
    for iSeq,length in enumerate(lengths):
        seq = seqBase[:length]
        writer.record("%000i_%i" % (iSeq,length),seq)


workDir = pjoin(options.testRunDir,"test_ImmAppInner")
makedir(workDir)
submitDir = pjoin(workDir,"submit")
makedir(submitDir)
seqDbPath = pjoin(workDir,"seqdb")
#seqDbLengths = n.logspace(n.log10(50000),n.log10(3000000),50,base=10).astype(int)
#seqDbLengths = n.linspace(50000,1000000,50).astype(int)
#seqDbLengths = [50000]*50 + [1000000]*50
seqDbLengths = [500000]*50
#genRandomSeqDbCat(seqDbPath,50000,20)
queryPath = pjoin(workDir,"query.fna")
#queryLengths = n.ones(10000,dtype=int)*5000
queryLengths = n.arange(5000,15001,1)

def genInput():
    genRandomSeqDbIid(seqDbPath,seqDbLengths)
    genRandomQueryZeroMarkov(queryPath,queryLengths)


if len(sys.argv) == 1:
    genInput()

opt = Struct()
opt.runMode = "batchDep"

# training opts

opt.seqDb = seqDbPath
seqDb = SeqDbFasta(opt.seqDb)
ids = seqDb.getIdList()
immIdToSeqIds = dict(((id,[id]) for id in ids))
immIdsFile = pjoin(workDir,"test.immapp.seqids.pkl")
dumpObj(immIdToSeqIds,immIdsFile)
opt.immIdToSeqIds = immIdsFile
opt.immDb = pjoin(workDir,"icm")

# scoring opts

opt.immIds = immIdsFile
opt.nImmBatches = 10
opt.inpSeq = queryPath
opt.outScoreComb = pjoin(workDir,"score.comb.pkl")
opt.outDir = pjoin(workDir,"score")
makedir(opt.outDir)


def train():
    os.chdir(submitDir)
    opt.mode = "train"
    imm = ImmApp(opt=opt)
    return imm.run()

def score():
    os.chdir(submitDir)
    opt.mode = "score"
    imm = ImmApp(opt=opt)
    return imm.run(depend=jobs)


if len(sys.argv) == 1:


    jobs = []

    jobs = train()
    jobs = score()

else:
    x = loadObj(opt.outScoreComb)
    for y in n.histogram(x.scores.argmax(1))[0]: print '*'*y
    s = x.scores
    sn = (s - s.mean(0))/s.std(0)
    for y in n.histogram(sn.argmax(1))[0]: print '*'*y
    import matplotlib as mpl
    mpl.use('agg')
    import pylab as pl
    from matplotlib import mlab
    mn = s.mean(1)
    vr = s.var(1)
    firstI = 1000
    firstY = vr[0:firstI*2].mean()
    firstX = x.lenSamp[firstI]
    #allX = n.arange(firstX,max(x.lenSamp),firstX)
    allX = x.lenSamp
    allY = allX.astype(float)/firstX*firstY
    #allY = (n.arange(len(allX))+1)*firstY
    pl.axes()
    pl.plot(x.lenSamp,mn,label="Observed")
    pl.title("Score Mean(Sample Length)")
    pl.legend()
    pl.savefig("samp_len_mean.png", format='png')
    pl.axes().clear()
    #pl.plot(x.lenSamp,s)
    #pl.savefig("samp_len.png", format='png')
    pl.plot(x.lenSamp,vr,"g.",label="Observed")
    wind=400
    mavgVr = mlab.movavg(vr,wind+1)
    pl.plot(x.lenSamp[wind/2:-wind/2],mavgVr,"b-",label="Movavg(Observed,400)",linewidth=3)
    pl.plot(allX,allY,"r-",linewidth=2,label="Expected as sum of i.i.d.")
    pl.title("Score Variance(Sample Length)")
    pl.legend()
    pl.savefig("samp_len_var.png", format='png')
