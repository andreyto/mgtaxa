"""Application interface to creation of features based on k-mer frequences out of biological sequences."""

from MGT.Shogun.Util import *
from shogun.Features import *
from MGT.Svm import *
from MGT.App import *
import functools as fn

__all__ = ["SeqRandomApp"]


class RandomSeqGen:

    def __call__(self,inSeq):
        pass

class RandomSeqGenShuffle:
    """Shuffle native sequences"""
    def __call__(self,inSeq):
        return nrnd.permutation(n.fromstring(inSeq,dtype='S1'))

class RandomSeqGenGcFixed:
    """GC content for each sequence is taken from a distribution, each sequence is shuffled"""
    def __init__(self,paramDistr):
        if paramDistr == "uniform":
            self.paramDistr = fn.partial(nrnd.uniform,0.2,0.8)
        else:
            raise ValueError(paramDistr)

    def __call__(self,inSeq):
        gc_p = self.paramDistr()
        nseq = len(inSeq)
        ngc = int(nseq*gc_p/2)
        nat = int(nseq*(1-gc_p)/2)
        return nrnd.permutation(n.fromstring('AT'*nat + 'CG'*ngc,dtype='S1'))

class RandomSeqGenGcFixedRepl:
    """GC content for each sequence is taken from a distribution, 
    each sequence is build by sampling with replacement from a shorter sequence"""

    def __init__(self,paramDistr):
        if paramDistr == "uniform":
            self.paramDistr = fn.partial(nrnd.uniform,20,80)
        else:
            raise ValueError(paramDistr)

    def __call__(self,inSeq):
        gc_ratio = self.paramDistr()
        alph = nrnd.permutation(n.fromstring('AT'*int(100-gc_ratio) + 'CG'*int(gc_ratio),dtype='S1'))
        #rec["feature"]=nrnd.permutation(n.resize(alph,len(rec["feature"]))).tostring()
        return alph[nrnd.randint(len(alph),size=len(inSeq))]

class RandomSeqGenGcMnom:
    """GC content is taken from a distribution and used to generate probabilities
    for a 4-nomial distribution, which generates nucleotides"""

    def __init__(self,paramDistr):
        if paramDistr == "uniform":
            self.paramDistr = fn.partial(nrnd.uniform,0.2,0.8)
        else:
            raise ValueError(paramDistr)
        self.alph = n.fromstring("ATCG",dtype='S1')

    def __call__(self,inSeq):
        gc_p = self.paramDistr()
        p = [gc_p/2,gc_p/2,(1-gc_p)/2,(1-gc_p)/2]
        nseq = len(inSeq)
        # nseq samples, with 1 experiment in each returns nseq x 4 array 
        # with just one '1' in each row, 
        # then argmax along rows gives the winning nucleotide index
        return self.alph[nrnd.multinomial(1,p,nseq).argmax(1)]

class RandomSeqGenAMnom:
    """'A' content is taken from a distribution and used to generate probabilities
    for a 4-nomial distribution, which generates nucleotides"""

    def __init__(self,paramDistr):
        if paramDistr == "uniform":
            self.paramDistr = fn.partial(nrnd.uniform,0.2,0.8)
        else:
            raise ValueError(paramDistr)
        self.alph = n.fromstring("ATCG",dtype='S1')

    def __call__(self,inSeq):
        gc_p = self.paramDistr()
        a_p = nrnd.uniform(0.1,0.4)
        p = [a_p,(1-a_p)/3,(1-a_p)/3,(1-a_p)/3]
        nseq = len(inSeq)
        # nseq samples, with 1 experiment in each returns nseq x 4 array 
        # with just one '1' in each row, 
        # then argmax along rows gives the winning nucleotide index
        return self.alph[nrnd.multinomial(1,p,nseq).argmax(1)]

class SeqRandomApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=("default",),
            dest="mode",default="default"),
            make_option("-i", "--in-seq",
            action="store", type="string",dest="inSeq"),
            make_option("-o", "--out-seq",
            action="store", type="string",dest="outSeq"),
            make_option("-b", "--balance",
            action="store", type="int",dest="balance",default=-2),
            make_option("-t", "--method",
            action="store", type="choice",choices=("shuffle","genGcFixed","genGcFixedRepl","genGcMnom","genAMnom"),
            dest="method",default="shuffle"),
            make_option(None, "--param-distr",
            action="store", type="choice",choices=("uniform","cluster"),
            dest="paramDistr",default="uniform"),
            make_option("-a", "--alphabet",
            action="store", type="choice",choices=("dna","protein"),
            dest="alphabet",default="dna"),
        ]
        return Struct(usage = "Generate random sequencies with various distributions of nucleotide frequencies\n"+\
                "%prog [options]",option_list=option_list)

    def doWork(self,**kw):
        opt = self.opt
        print "App options:\n", opt
        loadSeqPreproc = loadSeqPreprocIdent
        data = loadSeqs(opt.inSeq,preProc=loadSeqPreproc)

        print "Loaded " + showSvmDataCounts(data)
        if opt.balance >= -1:
            data = balance(data,opt.balance)
            print "Balanced to " + showSvmDataCounts(data)

        ## we need to get rid of long runs of degenerates,
        ## otherwise methods that shuffle the original sequence
        ## will create sequences that have no higher order k-mers at all
        if opt.alphabet == 'dna':
            applyToFeatData(data,transDegen)
        
        if opt.method == "shuffle":
            rndSeqGen = RandomSeqGenShuffle()
        elif opt.method.startswith("gen"):
            assert opt.alphabet == 'dna',"opt.method=%s is not implemented for alphabet " % (opt.method,opt.alphabet)
            if opt.method == "genGcFixed":
                rndSeqGen = RandomSeqGenGcFixed(paramDistr=opt.paramDistr)
            elif opt.method == "genGcFixedRepl":
                rndSeqGen = RandomSeqGenGcFixedRepl(paramDistr=opt.paramDistr)
            elif opt.method == "genGcMnom":
                rndSeqGen = RandomSeqGenGcMnom(paramDistr=opt.paramDistr)
            elif opt.method == "genAMnom":
                rndSeqGen = RandomSeqGenAMnom(paramDistr=opt.paramDistr)
            else:
                raise ValueError(opt.method)
        
        for rec in data:
            rec["feature"]=rndSeqGen(rec["feature"]).tostring()

        saveSeqs(data,opt.outSeq)
        print "Wrote " + showSvmDataCounts(data)


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(SeqRandomApp)

