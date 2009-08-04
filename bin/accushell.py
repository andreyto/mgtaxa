### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.PredProcessor import *
from MGT.Sql import *
from MGT.Debug import *
from glob import glob

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-s", "--split",
        action="store", type="int",dest="split",default=1),
        make_option("-d", "--decision-range",
        action="store", type="string",dest="decisionRange",default="-100,-99,10"),
    ]
    parser = OptionParser(usage = "usage: %prog [options] perf_file [perf_file ...]",option_list=option_list)
    (options, args) = parser.parse_args()
    return options,args

def loadTest(fileName):
    return numpy.asarray([ int(l.split(None,1)[0]) for l in open(fileName,'r') ],dtype='i4')

def labelToId(labToId,lab):
    return labToId[lab]

def predToSql(db,name,test,pred):
    arr = numpy.rec.fromarrays((test,pred),names="test,pred")
    predName = name+'_pred'
    db.createTableFromArray(name=predName,arr=arr,withData=True,returnInserter=False)

def perfToSql(db,labToId,name,perf):
    arr = numpy.rec.fromarrays((labelToId(labToId,perf.speLab),perf.spe),names="id,spe")
    predName = name+'_spe'
    db.createTableFromArray(name=predName,arr=arr,withData=True,returnInserter=False)
    arr = numpy.rec.fromarrays((labelToId(labToId,perf.senLab),perf.sen),names="id,sen")
    predName = name+'_sen'
    db.createTableFromArray(name=predName,arr=arr,withData=True,returnInserter=False)




#predFiles = sorted(glob('*.pred'))
#predFiles = [ "test-d-b-s_3-c_0.0008.pred" ]
#predFiles = [ "test-d-b-s_3-c_100.pred" ]
#predFiles = [ "test-d-b-s_3-c_0.0001.pred" ]
#predFiles = [ "test-c_8_g_0.000488281.pred" ]
#minDval = 0.3

options,args = getProgOptions()

predFiles = args


decisionRange = numpy.arange(*[ float(x) for x in options.decisionRange.split(',')])

db = DbSqlMy() 
test = loadTest("test.svm")
labToId = loadObj("labelToId.pkl")

if options.split > 1:
    sampSplit = loadObj("../sampSplit.test.pkl")
else:
    sampSplit = numpy.ones(len(test),dtype='i2')

for predFile in predFiles:
    print "Start pred file: %s" % predFile
    try:
        predRes = loadPred(predFile,len(test))
        predictor = Predictor1(predRes=predRes)
        pred = predRes.pred
        for minDval in decisionRange:
        #for minDval in numpy.arange(0,2,0.1):
            pred = predictor.predict(minDval=minDval)
            if options.split > 1:
                splitMask = sampSplit[:,options.split-2]
            else:
                splitMask = sampSplit
            for splitId in xrange(1,options.split+1):
                testSplit = test[splitMask == splitId]
                predSplit = pred[splitMask == splitId]
                perf = perfMetrics(test=testSplit,pred=predSplit)
                print "predFile=%s minDval=%s nSplits=%s iSplit=%s" % (predFile,minDval,options.split,splitId)
                print perf
                print "[lab:id:spe...]: ", "  ".join( [ "%i:%i:%.2f" % p for p in zip(perf.speLab,labelToId(labToId,perf.speLab),perf.spe) ] )
                print
                print "[lab:id:sen...]: ", "  ".join( [ "%i:%i:%.2f" % p for p in zip(perf.senLab,labelToId(labToId,perf.senLab),perf.sen) ] )
                #testId = labelToId(labToId,testSplit)
                #predId = labelToId(labToId,predSplit)
                perfToSql(db=db,labToId=labToId,name="tmp_perf",perf=perf)
                outM = open(predFile+'_dv-%s_spl-%s.cm.csv' % (minDval,splitId),'w')
                confMatrCsv(outM,perf.cm.m,labToId)
                outM.close()
                #debug_here()
    except:
        print "Error processing predictions file %s" % predFile
        raise
    print "End pred file: %s" % predFile

