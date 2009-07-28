"""Random access database to sequence chunks constructed as reference to sequence from SeqDb.
"""

from MGT.Common import *
from MGT.Taxa import *
from MGT.SeqDb import HdfSeqInd

from MGT.Hdf import *

from itertools import izip

import random

class HdfSampInd(pt.IsDescription):
    """HDF data type to describe one sequence sample.
    A sequence sample starts at some position in one sequence, possibly spans several
    full consequtive sequences and can end in an arbitrary position of the last sequence.
    Although, the sample length size is typically constant, we define it inside this type
    to allow for variable size samples in a more general case. E.g. this will cover WGS contigs
    for which we want to predict taxonomy without shredding them into contant size chunks.
    The indSeq, begin and sampLen fields are sufficient to pull the sample out from sequence database.
    The nSeq field is an optimization: it adds one byte to the size of the data structure.
    If the sample spans less than 256 sequences (including sequences where it starts and where it ends),
    nSeq holds that number. Otherwise, nSeq is set to zero. Knowing the number of sequences allows
    for more simple and faster Python code that pulls the sequence data.
    The reason for having this datatype at all instead of just creating another SeqInd is because
    we need to insert spacers when we pull the sample spanning several sequences."""

    indSeq        = pt.Int64Col(pos=1) # index of record in seq.ind
    # if begin + sample size > sequence size, sample continues in the next sequence(s)
    begin     = pt.Int64Col(pos=2) # offset of this sample's first element relative to the start of sequence
    nSeq      = pt.UInt8Col(pos=3) # number of sequences
    #not used yet:
    #sampLen   = pt.UInt32Col(pos=4) # length of sample

# I could not find in Numpy any methods to obtain numeric limits for integer datatypes (found only for floats).
# This converts -1 to the type of nSeq, which should give the max possible value for any unsigned integer type.
maxHdfSampInd_nSeq = HdfSampInd.columns['nSeq'].dtype.type(-1)
assert maxHdfSampInd_nSeq > 0, "HdfSampInd.nSeq must be unsigned type"
# Make it open range and make it work even for maximum machine unsigned int
maxHdfSampInd_nSeq = long(maxHdfSampInd_nSeq) + 1
#DEBUG:
#print "maxHdfSampInd_nSeq =", maxHdfSampInd_nSeq

class HdfTaxaSampInd(pt.IsDescription):
    id        = pt.Int32Col(pos=1) # taxonomy id
    # [begin,begin + size) form right-open-ended range, same as python range() or STL begin(),end()
    begin     = pt.Int64Col(pos=2) # index of first element within SampInd
    size       = pt.Int64Col(pos=3) # number of elements within SampInd that belong to this taxid


def hdfMakeSampIndexConcat(hdfFile,hdfGroup,hdfSeqInd,sampLen):
    """Create sample index that chunks sequence concatenated by taxid.
    @param hdfFile - create in this open HDFFile object
    @param hdfGroup - create in this group (group will be created if necessary)
    @param hdfSeqInd - reference sequence index data from this Table object.
    Index must list sequence in taxid order, and 'id' field set to taxid.
    @param sampLen - length of each sample
    @post hdfGroup has two tables created in it:
    - indSamp ( 
        indSeq - index of hdfSeqInd record where sample starts,
        begin - offset in sequence pointed by hdfSeqInd where sample starts
        )
    - indTaxa (
        id - taxonomy id
        begin - index of the first sample in indSamp with that taxonomy id
        size - number of consecutive records in indSamp with that taxonomy id
        )
    """
    lastSeqIndRec = hdfSeqInd.read(hdfSeqInd.nrows-1)[0]
    expectedSamples = long(lastSeqIndRec['begin']/sampLen)
    # This is almost certainly a gross overestimate:
    expectedTaxa = lastSeqIndRec['id']
    indSamp = hdfFile.createTable(hdfGroup, 
            'indSamp', 
            HdfSampInd, 
            "Sample Index",
            expectedrows=expectedSamples,
            createparents=True)
    indTaxa = hdfFile.createTable(hdfGroup, 
            'indTaxa', 
            HdfTaxaSampInd, 
            "Taxa Index",
            expectedrows=expectedTaxa,
            createparents=True)
    #size of internal memory buffer of this leaf
    indSamp.nrowsinbuf = 1024**2 #1M
    indTaxa.nrowsinbuf = 1024**2 #1M
    iterIndSeq = hdfSeqInd.iterrows()
    seqRec = iterIndSeq.next()
    #DEBUG:
    #print "seq",seqRec.nrow,seqRec
    sampRec = indSamp.row
    taxaRec = indTaxa.row
    sampBegin = 0L
    nSamp = 0L #somehow table.nrows is always zero till flush()
    taxaRec['id'] = seqRec['id']
    taxaRec['begin'] = 0
    taxaRec['size'] = 0
    def appendTaxaRec():
        taxaRec['size'] = nSamp - taxaRec['begin']
        if taxaRec['size'] > 0:
            #DEBUG:
            #indTaxa.flush()
            #print "taxa", taxaRec[:]
            taxaRec.append()
        return
    def nextSeqRec():
        seqRec = iterIndSeq.next()
        #DEBUG:
        #print "seq",seqRec.nrow,seqRec
        if seqRec.nrow % 100000 == 0:
            print "Processed %s sequence records" % seqRec.nrow
        if seqRec['id'] != taxaRec['id']:
            appendTaxaRec()
            taxaRec['id'] = seqRec['id']
            taxaRec['begin'] = nSamp
            return False
        else:
            return True
    # This implements a kind of "jumping frogs" algorithm -
    # each next sequence is split into sample chunks until
    # they cross sequence end, after which the last sample
    # is filled with consequtive sequences until a sequence
    # contains that sample's end.
    # A new taxonomy id seen in input resets the sample offset.
    try:
        while True:
            while True:
                sampRec['indSeq'] = seqRec.nrow
                sampRec['begin'] = sampBegin
                sampRec['nSeq'] = 1
                sampBegin += sampLen
                if sampBegin > seqRec['size']:
                    break
                #DEBUG:
                #indSamp.flush()
                #print "samp",indSamp.nrows, sampRec[:]
                sampRec.append()
                nSamp += 1
                if sampBegin == seqRec['size']:
                    nextSeqRec()
                    sampBegin = 0
            while True:
                sampBegin -= seqRec['size']
                if not nextSeqRec():
                    sampBegin = 0
                    break
                if sampBegin <= seqRec['size']:
                    nSeq = seqRec.nrow - sampRec['indSeq'] + 1
                    if nSeq >= maxHdfSampInd_nSeq:
                        nSeq = 0
                    sampRec['nSeq'] = nSeq
                    #DEBUG:
                    #indSamp.flush()
                    #print "samp",indSamp.nrows,sampRec[:]
                    sampRec.append()
                    nSamp += 1
                    if sampBegin == seqRec['size']:
                        nextSeqRec()
                        sampBegin = 0
                    break
    except StopIteration:
        appendTaxaRec()
    indSamp.flush()
    indTaxa.flush()


class HdfSampleReader(MGTOptions):

    def __init__(self,hdfSampFile,sampLen,spacer='N',featType="chararray"):
        """@param featType type for output feature, "array" - numpy.chararray('S1'), "string" - Python string"""
        MGTOptions.__init__(self)
        self.hdfSampFile = hdfSampFile
        self.sampLen = sampLen
        self.spacer = numpy.fromstring(spacer,'S1')
        self.featType = featType
        self._openHdf()
        taxaInd = {}
        for row in self.hdfTaxaInd:
            taxaInd[row['id']] = (row['begin'],row['size'])
        self.taxaInd = taxaInd
        self.clearSubSampler()

    def _openHdf(self):
        self.hdfFileSeq = pt.openFile(self.hdfSeqFile,mode="r")
        self.hdfSeq = self.hdfFileSeq.getNode(self.hdfSeqGroup,'seq')
        self.hdfFileActSeq = pt.openFile(self.hdfActSeqFile,mode="r")
        self.hdfSeqInd = self.hdfFileActSeq.getNode(self.hdfActSeqInd)
        self.hdfFileSamp = pt.openFile(self.hdfSampFile,mode="r")
        self.hdfSampInd = self.hdfFileSamp.getNode(self.hdfSampGroup,'indSamp')
        self.hdfTaxaInd = self.hdfFileSamp.getNode(self.hdfSampGroup,'indTaxa')

    def _closeHdf(self):
        self.hdfTaxaInd = None
        self.hdfSampInd = None
        self.hdfSeqInd = None
        self.hdfSeq = None
        if self.hdfFileSamp is not None:
            self.hdfFileSamp.close()
            self.hdfFileSamp = None
        if self.hdfFileActSeq is not None:
            self.hdfFileActSeq.close()
            self.hdfFileActSeq = None
        if self.hdfFileSeq is not None:
            self.hdfFileSeq.close()
            self.hdfFileSeq = None

    def getSampleSeq(self,sampInd):
        
        hdfSeqInd = self.hdfSeqInd
        hdfSeq = self.hdfSeq
        sampLen = self.sampLen
        spacer = self.spacer
        lenSpacer = len(spacer)
        nSeq = sampInd['nSeq']
        #print "DEBUG: 0: ", sampInd
        if nSeq == 1:
            seqInd = hdfSeqInd.read(start=sampInd['indSeq'],stop=sampInd['indSeq']+1)[0]
            #print "DEBUG: 00: ", seqInd
            retVal = hdfSeq.read(start=seqInd['begin']+sampInd['begin'],stop=seqInd['begin']+sampInd['begin']+sampLen)
        elif nSeq > 1:
            seqInds = hdfSeqInd.read(start=sampInd['indSeq'],stop=sampInd['indSeq']+nSeq)
            #print "DEBUG: 10: ", seqInds
            seqInds['begin'][0] += sampInd['begin']
            seqInds['size'][0] -= sampInd['begin']
            seqInds = numpy.column_stack((seqInds['begin'],seqInds['begin']+seqInds['size'],seqInds['size']))
            #print "DEBUG: 11: ", seqInds
            seqInds[-1,1] -= seqInds[:,2].sum() - sampLen
            #print "DEBUG: 12: ", seqInds
            chunks = [ hdfSeq.read(start=seqInds[i,0],stop=seqInds[i,1]) for i in xrange(nSeq) ]
            for i in xrange(1,2*len(chunks)-1,2):
                chunks.insert(i,spacer)
            retVal = numpy.concatenate(chunks)
        else:
            # We guess about the max array size for samp, because we do not know in advance
            # the number of spacers
            samp = numpy.zeros(sampLen*2,dtype='S1')
            # total accumulated sample length w/o spacers
            nSamp = 0
            # nSamp + all spacers added so far
            nSampTotal = 0
            sampSeqBegin = sampInd['begin']
            for seqInd in hdfSeqInd.iterrows(start=sampInd['indSeq'],stop=hdfSeqInd.nrows):
                #print "DEBUG: 20: ", seqInd
                sampSeqEnd = sampSeqBegin + sampLen - nSamp
                piece = hdfSeq.read(start=seqInd['begin'] + sampSeqBegin,stop=seqInd['begin']+min(sampSeqEnd,seqInd['size']))
                #print "DEBUG: 21: ", len(piece)
                samp[nSampTotal:nSampTotal+len(piece)] = piece
                nSampTotal += len(piece)
                if sampSeqEnd <= seqInd['size']:
                    break
                nSamp += len(piece)
                # add spacer
                samp[nSampTotal:nSampTotal+lenSpacer] = spacer
                nSampTotal += lenSpacer
                sampSeqBegin = 0
            retVal = samp[:nSampTotal]
        return retVal.view(n.chararray)


    def randomSamples(self,taxid,nSamples):
        (beginSamp,sizeSamp) = self.taxaInd[taxid]
        indSampSel = random.sample(xrange(beginSamp,beginSamp+sizeSamp),nSamples)
        for sampInd in self.hdfSampInd.itersequence(indSampSel):
            s = self.subSampler(self.getSampleSeq(sampInd))
            if self.featType == "string":
                s = s.tostring()
            yield dict(id=sampInd.nrow,feature=s)

    def setSubSampler(self,subSampler):
        self.subSampler = subSampler

    def clearSubSampler(self):
        self.subSampler = lambda sample: sample

    def setSubSamplerUniRandomEnd(self,minLen,maxLen):
        assert self.sampLen >= minLen and self.sampLen >= maxLen and minLen <= maxLen
        self.setSubSampler(SubSamplerUniRandomEnd(minLen=minLen,maxLen=maxLen))

    def checkFakeDb(self):
        seqChecker = getFakeSequenceChecker()
        seqChecker.loadGiTaxa()
        taxaInd = self.taxaInd
        taxind = taxaInd.items()
        taxind.sort(key=lambda item: item[1][0])
        hdfSampInd = self.hdfSampInd
        hdfSeqInd = self.hdfSeqInd
        hdfSeq = self.hdfSeq
        sampLen = self.sampLen
        iTaxa = 0
        iMismatch = 0
        for item in taxind:
            taxid = item[0]
            (beginSamp,sizeSamp) = item[1]
            for sampInd in hdfSampInd.itersequence(range(beginSamp,beginSamp+sizeSamp), sort=True):
                seqInd = hdfSeqInd.read(start=sampInd['indSeq'],stop=sampInd['indSeq']+1)[0]
                if taxid != seqInd['id']:
                    res = 'OK'
                    if taxid != seqInd['id']:
                        res = "!!"
                    print "%s iTaxa = %s taxaInd = %s sampInd = %s seqInd = %s taxid(seq[:40]) = %s" % \
                            (res,iTaxa,item,sampInd,seqInd,
                                    seqChecker.gi2taxa[int(hdfSeq.read(start=seqInd['begin'],stop=seqInd['begin']+40)\
                                            .tostring().split('x')[1])])
                    iMismatch += 1
            if iTaxa % 1000 == 0:
                print "Done %s taxa out of %s, with %s mismatches" % ( iTaxa, len(taxind), iMismatch )
            iTaxa += 1

        iTaxa = 0
        for item in taxind:
            taxid = item[0]
            nSamples = max(min(600,item[1][1]-2),1)
            for samp in self.randomSamples(taxid=taxid,nSamples=nSamples):
                seqChecker.checkArrayTaxid(taxid=taxid,seq=samp,spacer=self.spacer)
            if iTaxa % 1000 == 0:
                print "Done %s taxa out of %s" % ( iTaxa, len(taxind) )
            iTaxa += 1


    def checkDb(self):
        """Perform some HDF database consistency checks to make sure that we correctly created the data."""
        taxaInd = self.taxaInd
        hdfSampInd = self.hdfSampInd
        hdfSeqInd = self.hdfSeqInd
        hdfSeq = self.hdfSeq
        sampLen = self.sampLen
        taxaIndNew = {}
        hdfSampInd.nrowsinbuf = 1024**2
        hdfSeqInd.nrowsinbuf = 1024**2
        for sampInd in hdfSampInd.iterrows():
            # we just ignore samples that span too many sequences in this test
            assert sampInd['nSeq'] > 0
            seqInds = hdfSeqInd.read(start=sampInd['indSeq'],stop=sampInd['indSeq']+sampInd['nSeq'])
            taxidSamp = seqInds[0]['id']
            assert seqInds['size'].sum() - sampInd['begin'] >= sampLen
            assert (seqInds['id'] == taxidSamp).all()
            if taxidSamp in taxaIndNew:
                txind = taxaIndNew[taxidSamp]
                txind = (txind[0],txind[1]+1)
            else:
                txind = (sampInd.nrow,1)
            taxaIndNew[taxidSamp] = txind
            if sampInd.nrow % 100000 == 0:
                print "Processed %s samples out of %s."  % (sampInd.nrow,hdfSampInd.nrows) 

        assert taxaIndNew == taxaInd


    def close(self):
        self.taxaInd = None
        self._closeHdf()

    def getTaxaInd(self):
        return self.taxaInd

