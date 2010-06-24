### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Implementation of various classification tasks using SHOGUN library.
It relies on some patches to the SHOGUN itself, such as Platt's probability estimates."""

import numpy
from numpy import random as nrnd
import pdb

#from shogun.Kernel import * #GaussianKernel, WeightedDegreeStringKernel, WeightedCommWordStringKernel
from shogun.Features import Labels, SparseRealFeatures
from shogun.Classifier import PlattProb
#from shogun.PreProc import SortWordString, SortUlongString

def makeRandomSplits(nSamp,nSplits):
    """Generate a permuted index that splits nSamp samples into nSplits parts.
    @return array(nSamp) where each element is the split's index in the range [0,nSplits)"""
    n_ind = int(numpy.ceil(float(nSamp)/nSplits))
    splitIndex = nrnd.permutation(numpy.concatenate([numpy.arange(nSplits,dtype=int)]*n_ind))
    return splitIndex[:nSamp]

def splitTwo(features,labels,ratio,fileName):
    nSplits = int(1./ratio)
    nFeat = features.get_num_vectors()
    labVal = labels.get_labels()
    assert nSplits > 1
    splitIndex = makeRandomSplits(nFeat,nSplits)
    assert len(splitIndex) == nFeat
    splitCounts = numpy.bincount(splitIndex)
    assert sum(splitCounts>0) >= 2, "At least two non-empty splits are required"
    iSplit = 0
    assert splitCounts[iSplit] > 0
    indTs = numpy.where(splitIndex == iSplit)[0].astype('i4')
    indTr = numpy.where(splitIndex != iSplit)[0].astype('i4')
    lbTs = Labels(labVal[indTs])
    lbTr = Labels(labVal[indTr])
    ftTs = features.subsample(indTs)
    ftTr = features.subsample(indTr)
    print "Writing split with %d samples" % ftTr.get_num_vectors()
    ftTr.write_svmlight_file(fileName+'.1',lbTr)
    print "Writing split with %d samples" % ftTs.get_num_vectors()
    ftTs.write_svmlight_file(fileName+'.2',lbTs)

def crossValidate(svm,splitIndex=None,nSplits=None):
    """Perform a cross-validation.
    @param svm SVM object with labels and kernel assigned, which in turn has lhs features assigned
    @param splitIndex Numpy array with a split index for each corresponding feature vector
    @param nSplits If splitIndex is None, then feature vectors will be split randomly into nSplits
    @return Labels object for svm.kernel.lhs features assigned by cross-validation
    """
    assert splitIndex is None or ( nSplits is None and hasattr(splitIndex,"__len__") )
    kernel = svm.get_kernel()
    features = kernel.get_lhs()
    featuresRhs = kernel.get_rhs()
    labels = svm.get_labels()
    # numpy float array of training labels:
    labValTr = svm.get_labels().get_labels()
    labValTs = numpy.zeros_like(labValTr)
    nFeat = features.get_num_vectors()
    if nSplits is not None:
        splitIndex = makeRandomSplits(nFeat,nSplits)
    assert len(splitIndex) == nFeat
    splitCounts = numpy.bincount(splitIndex)
    assert sum(splitCounts>0) >= 2, "At least two non-empty splits are required for cross-validation"
    for iSplit in xrange(len(splitCounts)):
        if splitCounts[iSplit] > 0:
            indTs = numpy.where(splitIndex == iSplit)[0].astype('i4')
            indTr = numpy.where(splitIndex != iSplit)[0].astype('i4')
            ftTs = features.subsample(indTs)
            ftTr = features.subsample(indTr)
            lbVTr = labValTr[indTr]
            #print numpy.bincount(lbVTr.astype(int)+1)
            #pdb.set_trace()
            lbTr = Labels(lbVTr)
            kernel.init(ftTr,ftTr)
            svm.set_labels(lbTr)
            svm.set_batch_computation_enabled(True)
            svm.set_linadd_enabled(True)
            print "Training cross-validation split # %s" % iSplit
            svm.train()
            print 'Objective: %f num_sv: %d' % (svm.get_objective(), svm.get_num_support_vectors())
            kernel.init(ftTr, ftTs)
            svm.set_batch_computation_enabled(True)
            svm.set_linadd_enabled(True)
            print "Testing cross-validation split # %s" % iSplit
            labValTs[indTs] = svm.classify().get_labels()
    # Restore original state of kernel and svm objects
    kernel.init(features,featuresRhs)
    svm.set_labels(labels)
    return Labels(labValTs)


def makePlattProb(svm,splitIndex=None,nSplits=None):
    labCr = crossValidate(svm=svm,splitIndex=splitIndex,nSplits=nSplits)
    labTr = svm.get_labels()
    return PlattProb(labTr,labCr)

def convSparseToShog(data,delFeature=False):
    resFeat = SparseRealFeatures()
    resFeat.create_sparse_feature_matrix(len(data))
    for iRec in xrange(len(data)):
        feat = data[iRec]["feature"]
        resFeat.set_sparse_feature_vector(iRec,feat["ind"].astype('i4')-1,feat["val"].astype('f8'))
        if delFeature:
            data[iRec]["feature"] = None
    return resFeat

def selShogById(idLabs,data,featShog):
    """Subsample both data and corresponding Shogun feature object with IdLabels object.
    @param idLabs IdLabels instance
    @param data sample data record array
    @param featShog Shogun feature object with 1-to-1 correspondence to data, e.g. from convSparseToShog(data)
    @ret Struct(ind,data,featShog) where ind = idLabs.selDataInd(data), data is new data with "label" field updated
    from idLabs, featShog is new Shogun features"""
    assert len(data) == featShog.get_num_vectors()
    ind = idLabs.selDataInd(data)
    newData = data[ind]
    idLabs.setDataLab(newData)
    newFeatShog = featShog.subsample(ind.astype('i4'))
    return Struct(ind=ind,data=newData,featShog=newFeatShog)

