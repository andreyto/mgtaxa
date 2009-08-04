### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.SeqDb import hdfCheckActiveSeqInd
from MGT.SampDb import *

import time, random

printTiming = False
printData = True

hdfSampFile = "samp_1k/samp.hdf"
sampLen = 1000

reader = HdfSampleReader(hdfSampFile=hdfSampFile,sampLen=sampLen,spacer='    ')
taxaInd = reader.getTaxaInd()
#for taxid in taxaInd:
#    print taxid,taxaInd[taxid]


def testLocalSamples(taxid):
    nLoops = 1000
    nSamples = 3
    nSamplesTotal = nLoops*nSamples
    print "Starting repeated extraction of samples from one (short) taxa."
    start = time.time()
    for i in xrange(nLoops):
        for samp in reader.randomSamples(taxid=taxid,nSamples=nSamples):
            if printData:
                s = samp.tostring()
                print len(s), '"%s"' % s
    finish = time.time()
    timePerSamp = (finish-start)*1000/float(nSamplesTotal)
    if printTiming:
        print "Repeated one-taxa extraction: Done %s loops of %s samples, total of %s samples, %.3f ms per sample." % \
                (nLoops,nSamples,nSamplesTotal,timePerSamp)

def testRandomSamples():
    taxids = taxaInd.keys()
    #nrnd.shuffle(taxids)
    nLoops = 1000
    nSamples = 6
    nSamplesTotal = 0
    print "Starting repeated extraction of random samples from random taxa."
    start = time.time()
    for i in xrange(nLoops):
        taxid=taxids[random.randint(0,len(taxids)-1)]
        nSamplesActual = min(nSamples,taxaInd[taxid][1])
        nSamplesTotal += nSamplesActual
        for samp in reader.randomSamples(taxid=taxid,nSamples=nSamplesActual):
            if printData:
                s = samp.tostring()
                print len(s), '"%s"' % s
    finish = time.time()
    timePerSamp = (finish-start)*1000/float(nSamplesTotal)
    if printTiming:
        print "Repeated random taxa extraction: Done %s loops of max %s samples, total of %s samples, %.3f ms per sample." % \
                (nLoops,nSamples,nSamplesTotal,timePerSamp)

random.seed(1)
#testLocalSamples(taxid=262136)
#testRandomSamples()
## This will take several hours
#reader.checkDb()
#hdfCheckActiveSeqInd()
reader.checkFakeDb()
#reader.close()

