### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

"""Convert file in SVMLight format which has all features present into dense format"""

from MGT.PredProcessor import *
from MGT.Svm import SvmDenseFeatureWriterTxt

iArg=1

test = numpy.asarray([ int(l.split(None,1)[0]) for l in open(sys.argv[iArg],'r') ],dtype='i4')
iArg+=1

res = []
for predFile in sys.argv[iArg:]:
    svmOutFile = predFile+'.svm'
    print predFile, "->", svmOutFile
    try:
        pred,labels,indLab,dval = loadPred(predFile,len(test))
        svmWriter = SvmDenseFeatureWriterTxt(svmOutFile)
        for (label,values) in izip(test,dval):
            svmWriter.write(label=label,values=values)
            if svmWriter.numRec() % 10000 == 0:
                print "Wrote %s samples out of %s" % (svmWriter.numRec(),len(test))
    except Exception, msg:
        print "Error processing prediction file: ",
        print msg
        raise

