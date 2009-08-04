### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Kmers import *
from shogun.Features import SparseRealFeatures
import pdb

class KmerFeatures(KmerCounter):

    def __init__(self,kmerLen,revCompPolicy=RC_POLICY.MERGE):
        KmerCounter.__init__(self,kmerLen,revCompPolicy)

    def getSparseRealFeatures(self,sequences,method="frequences"):
        maxSeqLen = max( ( len(seq) for seq in sequences ) )
        kmer_ind = numpy.zeros(maxSeqLen,dtype='i8')
        if method == 'frequences':
            kmer_val = numpy.zeros(maxSeqLen,dtype='f4')
        else:
            kmer_val = numpy.zeros(maxSeqLen,dtype='i4')
        kmerMethod = getattr(self,method)
        resFeat = SparseRealFeatures()
        resFeat.create_sparse_feature_matrix(len(sequences))
        for iSeq in xrange(len(sequences)):
            seq = sequences[iSeq]
            if isinstance(seq,str):
                seq = numpy.fromstring(seq,'S1')
            self.process(seq)
            (size,total) = kmerMethod(kmer_val,kmer_ind)
            #print size, total, kmer_val[:10],kmer_ind[:10]
            resFeat.set_sparse_feature_vector(iSeq,kmer_ind[:size].astype('i4')-1,kmer_val[:size].astype('f8'))
        #pdb.set_trace()
        return resFeat

