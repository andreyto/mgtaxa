from MGT.SampDb import *
from MGT.Kmers import *

class KmerReader(HdfSampleReader,KmerCounter):

    def __init__(self,hdfSampFile,sampLen,kmerLen,spacer='N',maxDegen=None):
        if maxDegen is None:
            maxDegen = max(kmerLen*2,int(sampLen*0.05))
        self.maxDegen = maxDegen
        self.nDegenSkipped = 0
        HdfSampleReader.__init__(self,hdfSampFile=hdfSampFile,sampLen=sampLen,spacer=spacer)
        KmerCounter.__init__(self,kmerLen,RC_POLICY.MERGE)
        self.max_seqLen = sampLen * 2
        self.result = numpy.zeros(self.max_seqLen,dtype=numpy.dtype([('values','i4'),('indices','i8')]))
        self.resultF = numpy.zeros(self.max_seqLen,dtype=numpy.dtype([('values','f4'),('indices','i8')]))

    def _randomKmers(self,taxid,nSamples,method,result):
        self.nDegenSkipped = 0
        maxDegen = self.maxDegen
        for samp in self.randomSamples(taxid=taxid,nSamples=nSamples):
            self.process(samp)
            nDegen = self.sumDegenKmerCounts()
            (size,total) = method(result['values'],result['indices'])
            if nDegen < maxDegen:
                yield (result[:size],total)
            else:
                self.nDegenSkipped += 1
        
    def randomCounts(self,taxid,nSamples):
        for res in self._randomKmers(taxid,nSamples,self.counts,self.result):
            yield res

    def randomFrequences(self,taxid,nSamples):
        for res in self._randomKmers(taxid,nSamples,self.frequences,self.resultF):
            yield res

    def randomBits(self,taxid,nSamples):
        for res in self._randomKmers(taxid,nSamples,self.bits,self.result):
            yield res

    def randomSequences(self,taxid,nSamples):
        for res in self.randomSamples(taxid=taxid,nSamples=nSamples):
            yield (res,len(res))

    def randomFrequencesWriteSvmSparseTxt(self,label,taxidSamples,nSamples):
        maxDegen = self.maxDegen
        nWritten = 0
        for samp in self.randomSamples(taxid=taxidSamples,nSamples=nSamples):
            self.process(samp)
            #print "taxid = ", taxidSamples
            #print samp.tostring()
            #print "\n\n\n"
            #print "maxDegen = %s, sumDegenKmerCounts() = %s" % (maxDegen,self.sumDegenKmerCounts())
            nWritten += self.frequencesWriteSvmSparseTxt(label,maxDegen)
        return nWritten



class KmerReaderComb(HdfSampleReader):

    def __init__(self,hdfSampFile,sampLen,kmerLenRange,spacer='N',maxDegen=None):
        kmerLenRange = numpy.asarray(kmerLenRange)
        if maxDegen is None:
            maxDegen = max(kmerLenRange.max()*2,int(sampLen*0.05))
        self.kmerLenRange = kmerLenRange
        self.maxDegen = maxDegen
        self.nDegenSkipped = 0
        HdfSampleReader.__init__(self,hdfSampFile=hdfSampFile,sampLen=sampLen,spacer=spacer)
        self.kmerCounters = []
        firstIdState = 1
        # this will make each next kmer counter to assign
        # state IDs after the end of ID range used by the
        # previous counter
        for kmerLen in kmerLenRange:
            kc = KmerCounter(kmerLen,firstIdState)
            firstIdState = kc.getLastIdState()
            self.kmerCounters.append(kc)
        self.max_seqLen = sampLen * 2
        nKmerCounters = len(self.kmerCounters)
        self.result = numpy.zeros(self.max_seqLen*nKmerCounters,
                dtype=numpy.dtype([('values','i4'),('indices','i8')]))
        self.resultF = numpy.zeros(self.max_seqLen*nKmerCounters,
                dtype=numpy.dtype([('values','f4'),('indices','i8')]))

    def _randomKmers(self,taxid,nSamples,methodName,result):
        self.nDegenSkipped = 0
        maxDegen = self.maxDegen
        for samp in self.randomSamples(taxid=taxid,nSamples=nSamples):
            nDegen = 0
            size = 0
            total = 0
            for kmerCounter in self.kmerCounters:
                kmerCounter.process(samp)
                nDegen = max(nDegen,kmerCounter.sumDegenKmerCounts())
                method = getattr(kmerCounter,methodName)
                resCounter = result[size:]
                (sizeCounter,totalCounter) = method(resCounter['values'],resCounter['indices'])
                size += sizeCounter
                total += totalCounter
            if nDegen < maxDegen:
                yield (result[:size],total)
            else:
                self.nDegenSkipped += 1
        
    def randomCounts(self,taxid,nSamples):
        for res in self._randomKmers(taxid,nSamples,"counts",self.result):
            yield res

    def randomFrequences(self,taxid,nSamples):
        for res in self._randomKmers(taxid,nSamples,"frequences",self.resultF):
            yield res

    def randomBits(self,taxid,nSamples):
        for res in self._randomKmers(taxid,nSamples,"bits",self.result):
            yield res

