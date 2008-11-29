from MGT.Common import *

try:
    from MGTX.kmersx import *

    class KmerSparseFeatures(KmerCounter):

        def __init__(self,sampLen,kmerLen,rcPolicy=RC_POLICY.MERGE):
            KmerCounter.__init__(self,kmerLen,rcPolicy)
            self.max_seqLen = sampLen * 2
            self.result = numpy.zeros(self.max_seqLen,dtype=numpy.dtype([('values','i4'),('indices','i8')]))
            self.resultF = numpy.zeros(self.max_seqLen,dtype=numpy.dtype([('values','f4'),('indices','i8')]))

        
        def _kmers(self,samp,method,result):
            self.process(samp)
            (size,total) = method(result['values'],result['indices'])
            return result[:size]
            
        def kmerCounts(self,samp):
            return self._kmers(samp,self.counts,self.result)

        def kmerFrequences(self,samp):
            return self._kmers(samp,self.frequences,self.resultF)

        def kmerBits(self,samp):
            return self._kmers(samp,self.bits,self.result)

except ImportError:
    class SvmSparseFeatureWriterTxt:
        pass

