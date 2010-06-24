### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Common definitions used by feature-related modules"""
from MGT.Common import *

## Default Numpy dtype for label
labelDtype = n.dtype("f8")
## Default Numpy dtype for data split ids
splitDtype = n.dtype("i4")
## Default Numpy dtype for feature attribute index (used in sparse feature representation)
indDtype = n.dtype("i4")

## known formats of persistent storage of features
featIOFormats = ("txt","pkl")
## default value for format of persistent storage of features
defFeatIOFormat = "pkl"

class MGTSparseRealFeatures(object):
    """Datatype for sparse features of real numbers.
    Currently only used for Numpy datatype definitions,
    while the actual data is passed around as Numpy recarrays.
    @todo make it a subclass of numpy ndarray and hold the data here as well"""
    defDtype = n.dtype([("ind",indDtype),("val","f4")])

class MGTSparseIntFeatures(object):
    """Datatype for sparse features of integer numbers.
    Currently only used for Numpy datatype definitions,
    while the actual data is passed around as Numpy recarrays.
    @todo make it a subclass of numpy ndarray and hold the data here as well"""
    defDtype = n.dtype([("ind",indDtype),("val","i8")])

class MGTDenseRealFeatures(object):
    """Datatype for dense features of real numbers.
    Currently only used for Numpy datatype definitions,
    while the actual data is passed around as Numpy recarrays.
    @todo make it a subclass of numpy ndarray and hold the data here as well"""
    defDtype = n.dtype("f4")

class MGTSparseData(object):
    """Datatype for sparse data array that combines id,feature and label.
    Currently only used for Numpy datatype definitions,
    while the actual data is passed around as Numpy recarrays.
    @todo make it a subclass of numpy ndarray and hold the data here as well"""
    defDtype=[('label',labelDtype),('feature','O'),('id',idDtype)]
    
    @classmethod
    def makeEmpty(klass,numFeat):
        return n.empty(numFeat,dtype=klass.defDtype)

class RevCompl:
    def __init__(self):
        trans = ['N']*256
        for (c,o) in zip('ATCG','TAGC'):
            trans[ord(c)] = o
        for (c,o) in zip('atcg','tagc'):
            trans[ord(c)] = o
        self.trans = ''.join(trans)

    def __call__(self,s):
        return s.translate(self.trans)[::-1]

revCompl = RevCompl()

def featVecDenseToSparse(x):
    """Return sparse feature vector from a dense one"""
    ##@todo here we always save the index of last element, even if it is zero -
    ## otherwise when we convert later sparse feature into dense, we can lose
    ## dimensions at the end. This causes exceptions during prediction with 
    ## classifiers that operate on dense features.
    ## This has to be handled at some other level, but it is not clear where at this time.
    ind = x.nonzero()[0]
    if x[-1] == 0:
        ind = n.append(ind,len(x)-1)
    # sparse ind starts from 1
    return n.rec.fromarrays(((ind+1).astype(indDtype),x[(ind,)]),names=("ind","val"))

