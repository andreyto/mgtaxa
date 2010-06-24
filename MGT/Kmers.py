### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Svm import *
from MGTX.kmersx import *

## Feature normalization policy
class NORM_POLICY:
    ## no row-wise normalization
    NONE_ROW = 0x0001
    ## euclidian distance row-wise normalization (divide by sqrt(dot(x,x)))
    EU_ROW = 0x0002
    ## normalize by expected counts
    EXPECT = 0x0010
    ## normalize as frequency (n_k/sum(n_k) for each k-mer length). EU_ROW can be applied afterwards
    FREQ = 0x0020
    ## normalize by counts observed in a string-reversed sequence
    REVERSE = 0x0040
    

__all__ = ["NORM_POLICY","RC_POLICY","KmerSparseFeatures","KmerLadderSparseFeatures"]



class KmerSparseFeatures(KmerCounter):
    """Any of the kmerXXX methods may or may not return a reference to internal data,
    so make a copy if a guaranteed independent copy of the result is required."""

    def __init__(self,sampLen,kmerLen,rcPolicy=RC_POLICY.MERGE,normPolicy=None):
        if normPolicy is None:
            normPolicy = NORM_POLICY.EU_ROW
        KmerCounter.__init__(self,kmerLen,rcPolicy)
        self.normPolicy = normPolicy
        self.max_seqLen = sampLen * 2
        self.result = numpy.zeros(self.max_seqLen,dtype=MGTSparseIntFeatures.defDtype)
        self.resultF = numpy.zeros(self.max_seqLen,dtype=MGTSparseRealFeatures.defDtype)

    
    def _kmers(self,samp,method,result):
        if samp is not None:
            self.process(samp)
        (size,total) = method(result['val'],result['ind'])
        return result[:size]
        
    def kmerCounts(self,samp=None):
        """Return absolute k-mer counts.
        @param[in] samp If not None, process() will be called first for this data,
        and then count returned. Otherwise it will assume that process() has been called
        before (maybe multiple times), and simply finalize and return counts.
        """
        return self._kmers(samp,self.counts,self.result)

    def kmerFrequencies(self,samp=None):
        """Return frequencies normalized according to normalization policy supplied to ctor.
        @param[in] samp If not None, process() will be called first for this data,
        and then count returned. Otherwise it will assume that process() has been called
        before (maybe multiple times), and simply finalize and return counts.
        """
        res = self._kmers(samp,self.frequences,self.resultF)
        normPolicy = self.normPolicy
        assert normPolicy & NORM_POLICY.FREQ, normPolicy
        if normPolicy & NORM_POLICY.EU_ROW:
            res["val"] /= n.sqrt(n.dot(res["val"],res["val"]))
        elif normPolicy & NORM_POLICY.NONE_ROW:
            pass
        else:
            ValueError(normPolicy)
        return res


    def kmerBits(self,samp):
        return self._kmers(samp,self.bits,self.result)

class KmerLadderSparseFeatures(KmerCounterLadder):
    """Any of the kmerXXX methods may or may not return a reference to internal data,
    so make a copy if a guaranteed independent copy of the result is required."""

    def __init__(self,sampLen,kmerLen,rcPolicy=RC_POLICY.DIRECT,normPolicy=None):
        if normPolicy is None:
            normPolicy = NORM_POLICY.EXPECT | NORM_POLICY.NONE_ROW
        KmerCounterLadder.__init__(self,kmerLen,rcPolicy)
        self.normPolicy = normPolicy
        maxNumK = self.maxNumKmers(sampLen).sum()
        featDtype = MGTSparseRealFeatures.defDtype
        self.result = numpy.zeros(maxNumK,dtype=featDtype)
        self.valExp = numpy.zeros(maxNumK,dtype=featDtype["val"])
        self.sizes = numpy.zeros(kmerLen,dtype=featDtype["ind"])

    def kmerFrequencies(self,samp):
        if self.normPolicy & NORM_POLICY.REVERSE:
            self.normPolicy &= ~NORM_POLICY.REVERSE
            resObs = self.kmerFrequencies(samp).copy()
            valExp = self.valExp[:len(resObs)].copy()
            resRevObs = self.kmerFrequencies(samp[-1::-1])
            #valRevExp = self.valExp[:len(resRevObs)]
            self.normPolicy |= NORM_POLICY.REVERSE
            return sparseDiv(resObs,resRevObs,valExp)
        normPolicy = self.normPolicy
        self.process(samp)
        self.counts(self.result["val"],self.valExp,self.result["ind"],self.sizes)
        size = self.sizes.sum()
        res = self.result[:size]
        valE = self.valExp[:size]
        #TMP:
        # (X-E)/E
        # no need to normalize counts into frequency because
        # the frequency denominator N-k+1 will go away in the above equation
        assert not ((normPolicy & NORM_POLICY.EXPECT) and (normPolicy & NORM_POLICY.FREQ))
        if normPolicy & NORM_POLICY.EXPECT:
            res["val"] -= valE
            res["val"] /= valE
        elif normPolicy & NORM_POLICY.FREQ:
            vals = res["val"]
            iStart = 0
            for iEnd in self.sizes.cumsum():
                vals[iStart:iEnd] /= vals[iStart:iEnd].sum()
                iStart = iEnd
        else:
            ValueError(normPolicy)
        if normPolicy & NORM_POLICY.EU_ROW:
            res["val"] /= n.sqrt(n.dot(res["val"],res["val"]))
        elif normPolicy & NORM_POLICY.NONE_ROW:
            pass
        else:
            raise ValueError(normPolicy)
        return res

    def subFeatureSizes():
        return self.sizes


