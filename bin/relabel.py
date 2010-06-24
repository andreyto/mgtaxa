### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Relabel a set of SVMLight (training/testing) with consequitive labels.
The Numpy array that maps new labels to old labels is saved by pickle."""

import numpy
import sys
from cPickle import *

def loadLabels(inp,labSet):
    labSet.update([ int(l.split(None,1)[0]) for l in inp ])

iArgv = 1
indexOutFile = sys.argv[iArgv]
iArgv+=1

featureFiles = sys.argv[iArgv:]

labSet = set()

for ff in featureFiles:
    inp = open(ff,'r')
    loadLabels(inp,labSet)
    print "labSet: ", len(labSet), labSet
    inp.close()

labNew = dict(zip(sorted(labSet),range(1,len(labSet)+1)))

print "labNew: ", len(labNew), labNew

labToOld = numpy.zeros(len(labSet)+1,dtype='i4')

for (lOld,lNew) in labNew.iteritems():
    labToOld[lNew] = lOld

outNO = open(indexOutFile,'w')
dump(labToOld,outNO,-1)
outNO.close()

print "labToOld: ", len(labToOld), labToOld

for ff in featureFiles:
    inp = open(ff,'r')
    out = open(ff+'.new','w')
    for line in inp:
        lab, feature = line.split(None,1)
        lab = int(lab)
        labN = labNew[int(lab)]
        print "o: %i, n: %i" % (lab,labN)
        out.write("%i %s" % (labN,feature))
    out.close()

