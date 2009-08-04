### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.SampDbKmer import *

import time, random

printTiming = True
printData = False
getFrequences = True

hdfSampFile = "samp_1k/samp.hdf"
sampLen = 1000
kmerLen = 7

reader = KmerReader(hdfSampFile=hdfSampFile,sampLen=sampLen,kmerLen=kmerLen,spacer='    ')
taxaInd = reader.getTaxaInd()
#for taxid in taxaInd:
#    print taxid,taxaInd[taxid]

if getFrequences:
    randomMethod = reader.randomFrequences
else:
    randomMethod = reader.randomVectors

def testLocalSamples(taxid):
    nLoops = 1000
    nSamples = 3
    nSamplesTotal = nLoops*nSamples
    print "Starting repeated extraction of samples from one (short) taxa."
    start = time.time()
    for i in xrange(nLoops):
        for (samp,total) in randomMethod(taxid=taxid,nSamples=nSamples):
            if printData:
                s = samp
                print len(s), total, '"%s"' % s
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
        for samp,total in randomMethod(taxid=taxid,nSamples=nSamplesActual):
            if printData:
                s = samp
                print len(s), total, '"%s"' % s
    finish = time.time()
    timePerSamp = (finish-start)*1000/float(nSamplesTotal)
    if printTiming:
        print "Repeated random taxa extraction: Done %s loops of max %s samples, total of %s samples, %.3f ms per sample." % \
                (nLoops,nSamples,nSamplesTotal,timePerSamp)

random.seed(1)
testLocalSamples(taxid=262136)
testRandomSamples()
#reader.checkDb()
reader.close()

