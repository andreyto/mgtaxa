### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.PredProcessor import *

iArg = 1
if sys.argv[iArg] == "-s":
    doPredict = False
    iArg+=1
else:
    doPredict = True
    startDval = float(sys.argv[iArg])
    iArg+=1
    endDval = float(sys.argv[iArg])
    iArg+=1
    stepDval = float(sys.argv[iArg])
    iArg+=1

test = numpy.asarray([ int(l.split(None,1)[0]) for l in open(sys.argv[iArg],'r') ],dtype='i4')
iArg+=1

iPredClass = 1
predClass = Predictor1
if sys.argv[iArg] == "-p":
    iArg+=1
    iPredClass = int(sys.argv[iArg])
    iArg+=1
    print "Using alternative predictor class #%s" % iPredClass
    predClass = eval("Predictor%s" % iPredClass)

labelConv = None
rankConv = "class"
res = []
for predFile in sys.argv[iArg:]:
    print predFile
    try:
        resFile = []
        #pred = numpy.asarray([ int(l) for l in open(predFile,'r') ],dtype=int)
        predRes = loadPred(predFile,len(test))
        if False:
            if labelConv is None:
                labelConv = LabelConv(labels=predictor.labels,indLab=predictor.indLab)
                test = labelConv.convSampLabels(lab=test,level=rankConv)
            predictorL = PredictorL(labelConv=labelConv,predRes=predRes)
        if doPredict:
            predictor = predClass(predRes=predRes)
            for minDval in arange(startDval,endDval,stepDval):
                pred = predictor.predict(minDval=minDval)
                #pred = predictorL.predict(pred=pred,level=rankConv,minDval=minDval)
                perf = perfMetrics(test=test,pred=pred)
                resFile.append((minDval,around(perf.senMean*100),around(perf.speMean*100)))
        else:
            perf = perfMetrics(test=test,pred=predRes.pred)
            resFile.append((around(perf.senMean*100),around(perf.speMean*100)))
        res.append([predFile] + list(resFile))
    except Exception, msg:
        print "Error processing prediction file: ",
        print msg
        raise
print '\n'.join(["ROC: %s" % (row,) for row in res])

