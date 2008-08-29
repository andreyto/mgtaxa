"""Prepare a string of class weights for SVMLIB/LIBLINEAR training.
Quoting from LIBLINEAR Readme:
> train -c 10 -w1 2 -w2 5 -w3 2 four_class_data_file

Train four classifiers:
    positive        negative        Cp      Cn
    class 1         class 2,3,4.    20      10
    class 2         class 1,3,4.    50      10
    class 3         class 1,2,4.    20      10
    class 4         class 1,2,3.    10      10
"""

import numpy
import sys

def loadLabels(fileName):
    return numpy.asarray([ int(l.split(None,1)[0]) for l in open(fileName,'r') ])

iArgv = 1
svmTrFile = sys.argv[iArgv]
iArgv+=1

lab = loadLabels(svmTrFile)

labcnt = numpy.bincount(lab)

labcntMax = numpy.max(labcnt)

indW = numpy.where(labcnt > 0)

print [ c for c in labcnt[indW] ]
print indW[0]

weight = (float(labcntMax)/labcnt[indW])**2

#weight[weight>=10] = 10.

print ' '.join([ "-w%s %.3f" % tuple(pair) \
        for pair in numpy.rec.fromarrays((indW[0],weight),names='ind,weight')])


