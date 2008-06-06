from MGT.SampDb import *
from MGTX.kmersx import *

class KmerReader(HdfSampleReader,KmerCounter):

    def __init__(self,hdfSampFile,sampLen,kmerLen,spacer='N'):
        HdfSampleReader.__init__(self,hdfSampFile=hdfSampFile,sampLen=sampLen,spacer=spacer)
        KmerCounter.__init__(self,kmerLen,True)
        self.max_seqLen = sampLen * 2
        self.result = numpy.zeros(self.max_seqLen,dtype=numpy.dtype([('values','i4'),('indices','i8')]))
        self.resultF = numpy.zeros(self.max_seqLen,dtype=numpy.dtype([('values','f4'),('indices','i8')]))

    def randomCounts(self,taxid,nSamples):
        result = self.result
        for samp in self.randomSamples(taxid=taxid,nSamples=nSamples):
            self.process(samp)
            (size,total) = self.counts(result['values'],result['indices'])
            yield (result[:size],total)
        
    def randomFrequences(self,taxid,nSamples):
        result = self.resultF
        for samp in self.randomSamples(taxid=taxid,nSamples=nSamples):
            self.process(samp)
            (size,total) = self.frequences(result['values'],result['indices'])
            yield (result[:size],total)

    def randomFrequencesWriteSvmSparseTxt(self,taxid,nSamples):
        for samp in self.randomSamples(taxid=taxid,nSamples=nSamples):
            self.process(samp)
            self.frequencesWriteSvmSparseTxt(taxid)

