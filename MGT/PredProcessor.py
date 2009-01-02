from MGT.Util import *
from MGT.Taxa import *
import sys, numpy
from numpy import *
from itertools import izip
from cPickle import dump, load
import pdb
import os

def confusionMatrix(maxLabel,tests):
    m = numpy.zeros((maxLabel,maxLabel),'i4')
    for test in tests:
        m[test[0],test[1]] += 1
    # Specificity is sensitive to unbalanced
    # testing set (small classes will be swamped
    # by false positives from large classes).
    # We balance the matrix by multiplying
    # rows with weights that equalize the number
    # of test cases across classes (as though we
    # repeated the tests for small classes in
    # order to match the larger ones)
    nT = m.sum(-1)
    maxNT = numpy.max(nT)
    w = nT.astype('f4') / maxNT
    w[w == 0] = 1
    w = (1. / w)
    w[w == 0] = 1
    mb = numpy.around((m.transpose() * w).transpose()).astype('i4')
    return Struct(m=m,mb=mb)


class PerfMetrics(Struct):

    def setLabToName(self,labToName):
        self.labToName = labToName
    
    def toNameSpe(self):
        return [ (self.labToName[spe["label"]],spe) for spe in self.spe ]

    def toNameSpeStr(self):
        return "[name:lab:spe...]: ", "  ".join( [ "%s:%s:%.2f" % (p[0],p[1]["label"],p[1]["val"]) for p in  self.toNameSpe() ] )

    def toNameSen(self):
        return [ (self.labToName[sen["label"]],sen) for sen in self.sen ]

    def toNameSenStr(self):
        return "[name:lab:sen...]: ", "  ".join( [ "%s:%s:%.2f" % (p[0],p[1]["label"],p[1]["val"]) for p in  self.toNameSen() ] )

    def confMatrCsv(self,out,m=None):
        if m is None:
            m = self.cm.m
        if isinstance(out,str):
            out = open(out,'w')
            closeOut = True
        else:
            closeOut = False
        labToId = self.labToName
        lab = numpy.arange(len(m),dtype='i4')
        name = [ labToName[l] for l in lab ]
        out.write('0,0,')
        for l in lab: out.write("%i," % l)
        out.write('\n')
        out.write('0,0,')
        for i in name: out.write("%i," % i)
        out.write('\n')
        for l in lab:
            out.write("%i,%i," % (l,name[l]))
            for x in m[l]: out.write("%i," % x)
            out.write('\n')
        if closeOut:
            out.close()

def perfMetrics(test,pred,balanceCounts=True):
    maxLabel = max(test.max(),pred.max()) + 1
    cm = confusionMatrix(maxLabel,izip(test,pred))
    if balanceCounts:
        m = cm.mb
    else:
        m = cm.m
    #print m
    #print cm.m
    #print cm.mb
    mr = m[1:,1:]
    t = m[1:,:].sum(axis=1)
    p = mr.sum(axis=0)
    tp = mr.diagonal()
    tNzW = numpy.where(t > 0)
    senT = t[tNzW]
    sen = (tp[tNzW].astype('f4')/t[tNzW])
    senLab = tNzW[0] + 1 #because mr = m[1:,1:]
    sen = n.rec.fromarrays([senLab,sen],names="label,val")
    senMean = sen["val"].mean()
    pNzW = numpy.where(p > 0)
    spe = (tp[pNzW].astype('f4')/p[pNzW])
    speLab = pNzW[0] + 1 #because mr = m[1:,1:]
    spe = n.rec.fromarrays([speLab,spe],names="label,val")
    if len(spe) > 0:
        speMean = spe["val"].mean()
        speMin = spe["val"].min()
        speMeanTP = spe["val"][spe["val"]>0].mean()
    else:
        speMean = 0
        speMin = 0
        speMeanTP = 0
    acc = float(tp.sum())/t.sum()
    #print "Spe: ", spe
    #print "Sen: ", sen
    print \
            "TP: "+`tp.sum()`+" T: "+`t.sum()`+" P: "+`p.sum()`+" S: "+`m.sum()`+\
            " CT: "+`(t>0).sum()`+" CP: "+`(p>0).sum()`+" CTP: "+`(tp>0).sum()`+\
            " Sen: %.2f"%(senMean*100)+" Spe: %.f"%(speMean*100)+" SpeMin: %.f"%(speMin*100)+\
            " SpeTP: %.f"%(speMeanTP*100)+" Acc: %.f"%(acc*100)
    pm = PerfMetrics(cm=cm,sen=sen,spe=spe,senMean=senMean,speMean=speMean,
            speMin=speMin,speMeanTP=speMeanTP,acc=acc)
    return pm

def loadPred(fileName,nTests):
    pklFileName = fileName+'.pkl'
    if os.path.isfile(pklFileName):
        inp = open(pklFileName,'r')
        res = load(inp)
        if res[2] is not None:
            res[2][0] = 999999999
        inp.close()
    else:
        inp = open(fileName,'r')
        line = inp.readline()
        labels = None
        if line.startswith("labels"):
            labels = numpy.fromstring(line[len('labels'):],dtype='i4',sep=' ')
            pred = numpy.zeros(nTests,dtype='i4')
            dval = numpy.zeros((nTests,len(labels)),dtype='f4')
            iPred = 0
            for line in inp:
                val = numpy.fromstring(line,dtype='f4',sep=' ')
                pred[iPred] = int(val[0])
                dval[iPred] = val[1:]
                iPred += 1
                if iPred % 10000 == 0:
                    print "Loaded %s/%s predictions" % (iPred,len(pred))
            assert iPred == nTests, "Number of tests = %s, number of predictions = %s" % (nTests,iPred)
            indLab = numpy.zeros(labels.max()+1,dtype='i4')
            indLab[:] = 999999999
            for i in xrange(len(labels)): indLab[labels[i]] = i
            res = (pred,labels,indLab,dval)
        else:
            inp.seek(0)
            pred = numpy.fromfile(inp,dtype='i4',sep='\n')
            iPred = len(pred)
            assert iPred == nTests, "Number of tests = %s, number of predictions = %s" % (nTests,iPred)
            res = (pred,None,None,None)
        inp.close()
        try:
            pout = open(pklFileName,'w')
            dump(res,pout,-1)
            pout.close()
        except MemoryError:
            try:
                pout.close()
            except:
                pass
            if os.path.isfile(pklFileName):
                os.remove(pklFileName)
    return Struct(pred=res[0],labels=res[1],indLab=res[2],dval=res[3])


class Predictor:
    
    def __init__(self,predRes):
        self.pred = predRes.pred
        self.labels = predRes.labels
        self.indLab = predRes.indLab
        self.dval = predRes.dval


class Predictor0(Predictor):
    """A do-nothing predictor"""
    
    def predict(self,minDval):
        return self.pred


class Predictor1(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)

    def predict(self,minDval):
        pred = self.pred
        dval = self.dval
        indLab = self.indLab
        dvalPred = numpy.asarray([dval[tuple(i)] \
                for i in izip(numpy.arange(len(dval),dtype='i4'),indLab[pred])],dtype='f4')
        predNew = (dvalPred>minDval).choose(0,pred)
        #predNew = pred
        return predNew

class Predictor2(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)
        #iDval = numpy.argsort(dval,-1)
        self.dval.sort(-1)

    def predict(self,minDval):
        dval = self.dval
        pred = self.pred
        predNew = numpy.logical_and(dval[:,-1] > 0,dval[:,-2] < 0).choose(0,pred)
        return predNew

class Predictor3(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)
        #iDval = numpy.argsort(dval,-1)
        self.dval.sort(-1)

    def predict(self,minDval):
        dval = self.dval
        pred = self.pred
        predNew = (dval[:,-1] - dval[:,-2] > minDval).choose(0,pred)
        return predNew

class Predictor4(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)
        #iDval = numpy.argsort(dval,-1)
        self.dval.sort(-1)

    def predict(self,minDval):
        dval = self.dval
        pred = self.pred
        predNew = (numpy.abs((dval[:,-1] - dval[:,-2])/(dval[:,-1] + dval[:,-2])) > minDval).choose(0,pred)
        return predNew


class LabelConv:
    """Convert labels into ranks in the lineage."""

    def __init__(self,labels,indLab):
        self.labels = labels
        self.indLab = indLab
        self.levels = TaxaLevels()
        print "Loading taxonomy tree"
        store = NodeStoragePickleSep(fileName="/home/atovtchi/work/mgtdata/samp_1k/taxaTree.postMarkTraining.pkl")
        self.taxaTree = TaxaTree(storage=store)
        print "Taxonomy tree loaded"
        self.levelNames = self.levels.getLevelNames()
        self.labelToId = loadObj("labelToId.pkl")
        self._indexDecisionColumnsToLevels(labels)

    def convSampLabels(self,lab,level):
        iLev = self.levelNames.index(level)
        print sorted(self.labLevName[iLev].items())
        return self.labToLevLab[lab,iLev]

    def getTaxaName(self,taxid):
        if taxid == 0:
            return "Reject"
        else:
            return self.taxaTree.getNode(taxid).name

    def _indexDecisionColumnsToLevels(self,labels):
        taxaTree = self.taxaTree
        colToLev = numpy.zeros((len(labels),len(self.levelNames)),dtype='i4')
        labelToId = self.labelToId
        levels = self.levels
        taxaTree = self.taxaTree
        indLab = self.indLab
        for iCol in xrange(len(labels)):
            taxid = labelToId[labels[iCol]]
            colToLev[iCol,:] = levels.lineageFixedList(node=taxaTree.getNode(taxid),null=0)
        colToLevLab = numpy.zeros_like(colToLev)
        # levLab[level ind] = { taxid -> label }
        self.levLab = []
        # labLev[level ind] = { label -> taxid }
        self.labLev = []
        # same but { label -> taxa name }
        self.labLevName = []
        for iLev in xrange(colToLev.shape[1]):
            uniqLev = numpy.unique(colToLev[:,iLev])
            if uniqLev[0] == 0:
                uniqLev = uniqLev[1:]
            levLab = dict(izip(uniqLev,numpy.arange(1,len(uniqLev)+1,dtype='i4')))
            levLab[0] = 0
            for iCol in xrange(len(colToLev)):
                colToLevLab[iCol,iLev] = levLab[colToLev[iCol,iLev]]
            self.levLab.append(levLab)
            self.labLev.append(dict([ (item[1],item[0]) for item in levLab.items() ]))
            self.labLevName.append(dict([ (item[1],self.getTaxaName(item[0])) 
                for item in levLab.items() ]))
        sh = colToLevLab.shape
        labToLevLab = numpy.zeros((sh[0]+1,sh[1]),dtype='i4')
        labToLevLab[labels,:] = colToLevLab[indLab[labels],:]
        self.colToLev = colToLev
        self.colToLevLab = colToLevLab
        self.labToLevLab = labToLevLab


class PredictorL(Predictor):

    def __init__(self,labelConv,**kw):
        Predictor.__init__(self,**kw)
        self.labelConv = labelConv
        self.levelNames = labelConv.levelNames
        self.colToLev = labelConv.colToLev

    def predict(self,pred,minDval,level):
        dval = self.dval
        indLab = self.indLab
        return self.labelConv.convSampLabels(lab=pred,level=level)

