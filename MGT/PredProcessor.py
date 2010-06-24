### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Taxa import *
from itertools import izip
import pdb
import os

def confusionMatrix(maxLabel,tests):
    """Build a confusion matrix from (true label, predicted label) pairs.
    @param maxLabel upper bound on the label value (max value + 1)
    @param tests sequence of pairs (true label, predicted label)
    @return {m=raw confusion matrix,mb=class balanced confusion matrix} where rows are positives and columns are predictions"""
    m = n.zeros((maxLabel,maxLabel),'i4')
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
    maxNT = n.max(nT)
    w = nT.astype('f4') / maxNT
    w[w == 0] = 1
    w = (1. / w)
    w[w == 0] = 1
    mb = n.around((m.transpose() * w).transpose()).astype('i4')
    return Struct(m=m,mb=mb)


class PerfMetrics(Struct):
    """Calculates performance metrics from pairs of true and predicted labels for each test case.
    Designated reject group is treated in a special way."""

    ## This label is designated for reject group.
    ## Reject group contributes to sensitivity but has undefined specificity
    rejectLabel = 0

    
    def __init__(self,test,pred,balanceCounts=True,maxLabel=None):
        """Constructor.
        @param test sequence of test (true) labels (int)
        @param pred sequence of predicted labels (int)
        @param balanceCounts if True, use confusion matrix balanced by per-class test counts
        @param maxLabel maximum value for labels plus one, if None - computed from actual test and pred values
        @post this object has masked array attributes for specificity, sensitivity, true positives etc for each class
        @todo Add additional metrics. One particular problem with per class specificities is when we
        have no testing samples for a given class, then even one false positive will create
        per class specificity of zero, weighing down the average specificity.
        @todo check out MaskedArrayAlternative in scipy, which has (experimental) masked recod array - 
        this way we can create a single record array (table) with fields like sen,spe,tp etc and one record
        for each class, and have each field in each record masked separately. Sort of like Excel table with
        possible empty values."""
        if maxLabel is None:
            maxLabel = max(test.max(),pred.max()) + 1
        cm = confusionMatrix(maxLabel,izip(test,pred))
        if balanceCounts:
            m = cm.mb
        else:
            m = cm.m
        #mr = m[1:,1:]
        #t = m[1:,:].sum(axis=1)
        mr = m[:,:]
        t = m[:,:].sum(axis=1)
        p = mr.sum(axis=0)
        tp = mr.diagonal()
        tMask = (t == 0)
        tm = nma.masked_array(t,mask=tMask)
        tpf = tp.astype('f4')
        sen = tpf/tm
        senMean = sen.mean()
        senMin = sen.min()
        senStd = sen.std()
        pMask = (p == 0)
        # For reject group (label 0) specificity is always undefined
        pMask[self.rejectLabel] = True
        pm = nma.masked_array(p,mask=pMask)
        spe = tpf/pm
        speMean = spe.mean()
        speMin = spe.min()
        speStd = spe.std()
        acc = tpf.sum()/t.sum()
        self.cm = cm
        self.sen = sen
        self.spe = spe
        self.senMean = senMean
        self.speMean = speMean
        self.senMin = senMin
        self.speMin = speMin
        self.senStd = senStd
        self.speStd = speStd
        self.acc    = acc
        self.tp = tpf
        self.p = pm
        self.t = tm
        self.names = n.arange(len(self.t))
        self.labToName = dict(enumerate(self.names))
        #if options.debug >= 1:
        #    self.summary()

    def setLabToName(self,labToName):
        self.labToName = labToName
        # create masked array out of dict
        self.names = dictToMaskedArray(labToName,maxInd=len(self.p),dtype='S32')

    def summary(self,out=None):
        if out is None:
            out = sys.stdout
        metrNames = ("names","spe","sen","tp","t","p")
        metrDict = dict([(metrName,getattr(self,metrName)) for metrName in metrNames])
        for name in ("spe","sen"):
            metrDict[name] = floatsToPcntStr(metrDict[name])
        metrVals = zip(*[metrDict[name] for name in metrNames])
        out.write(' '.join(["%10s" % x for x in metrNames])+'\n')
        out.write(' '.join(["%10s" % x for x in ["mean"]+floatsToPcntStr([self.speMean,self.senMean])])+'\n')
        out.write(' '.join(["%10s" % x for x in ["min"]+floatsToPcntStr([self.speMin,self.senMin])])+'\n')
        out.write(' '.join(["%10s" % x for x in ["std"]+floatsToPcntStr([self.speStd,self.senStd])])+'\n')
        for vals in metrVals:
            out.write(' '.join(["%10s" % x for x in vals])+'\n')


    def confMatrCsv(self,out,m=None):
        """Write confusion matrix in CSV format"""
        if m is None:
            m = self.cm.m
        if isinstance(out,str):
            out = open(out,'w')
            closeOut = True
        else:
            closeOut = False
        labToName = self.labToName
        lab = n.arange(len(m),dtype='i4')
        name = [ labToName.get(l,"") for l in lab ]
        out.write('0,0,')
        for l in lab: out.write("%i," % l)
        out.write('\n')
        out.write('0,0,')
        for s in name: out.write("%s," % s)
        out.write('\n')
        for l in lab:
            out.write("%i,%s," % (l,name[l]))
            for x in m[l]: out.write("%i," % x)
            out.write('\n')
        if closeOut:
            out.close()

    def getMetricsNames(self):
        return self.keys()



class PerfMetricsSet:
    """Holds multiple sets of performance metrics for the same samples and different sets of parameters.
    Can return selected metrics along with parameter values as numpy record arrays."""

    def __init__(self,param,perf):
        """Constructor.
        @param param recarray[N_param]
        @param perf sequence PerfMetrics[N_param]
        """
        self.param = param
        self.perf = perf

    def getMetrics(self,names):
        """Return a numpy recarray with dtype ["param",param_dtype),("val",[("name1","type1"),...])] and size N_param.
        @todo Currently all metrics types will be converted to 32 bit floats."""
        val = n.rec.fromrecords([ tuple([ getattr(pm,name) for name in names ]) for pm in self.perf ],
                names=','.join(names),formats=','.join(["f4"]*len(names)))
        if self.param is not None:
            fields = ("param","val")
            arrs = (self.param,val)
        else:
            fields = ("val",)
            arrs = (val,)
        return recFromArrays(arrs,names=fields)

    def exportMetricsCsv(self,names,out):
        """Call getMetrics(names) and export the result as CSV file"""
        arr = self.getMetrics(names=names)
        saveRecArrayAsCsv(arr,out=out,sep='\t')

    def getMetricsNames(self):
        return self.perf[0].getMetricsNames()

    def joinParam(self,param):
        """Join existing parameter records with the new one(s).
        That means append new fields to the parameter records.
        @param param a single recarray.record or recarray[N_param]
        """
        if self.param is not None:
            self.param = joinRecArrays([param,self.param])
        else:
            newpar = n.zeros(len(self.perf),dtype=param.dtype)
            newpar[:] = param
            self.param = newpar

    @classmethod
    def concatenate(klass,sets):
        """Return a union of PerfMetricsSets by concatenating their data.
        @param a sequence of PerfMetricsSet objects with identical param fields
        @return a union PerfMetricsSet object.
        This is a class method, so even if it is called as self.concatenate(sets),
        the content of 'self' is not used.
        """
        return klass(param=n.concatenate([s.param for s in sets]),
                perf=n.concatenate([s.perf for s in sets]))


class Predictions:
    """Holds multiple sets of predicted labels for the same samples and different sets of parameters.
    Also holds associated sample IDs and true labels, and can compute performance metrics."""

    def __init__(self,labPred,param,idPred):
        """Constructor.
        @param labPred array[N_param,N_samp]
        @param param recarray[N_param]
        @param idPred recarray("id","label")[N_samp] where "label" is true label (arbitrary value if not known)
        """
        self.labPred = labPred
        self.param = param
        self.idPred = idPred

    def calcPerfMetrics(self,idLab,confMatrFileStem=None,keepConfMatr=False,balanceCounts=False):
        """Compute and return performance metrics.
        @param idLab IdLabels (currently true label values are taken from self.idPred array rather than from idLab)
        @param confMatrFileStem if not None, save confusion matrices for each prediction set in files with this stem name
        @param keepConfMatr if True, keep the confusion matrix data inside the returned value, otherwise delete it
        """
        perf = []
        for iPred in range(len(self.labPred)):
            labP = self.labPred[iPred]
            pm = PerfMetrics(self.idPred["label"].astype(int),labP,balanceCounts=balanceCounts)
            pm.setLabToName(idLab.getLabToName())
            if confMatrFileStem is not None:
                pm.confMatrCsv(confMatrFileStem+'.%00i.cm.csv' % iPred)
            if not keepConfMatr:
                delattr(pm,"cm")
            perf.append(pm)
            if options.debug > 0:
                print self.param[iPred]
                pm.summary()
        return PerfMetricsSet(perf=perf,param=self.param)

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
        dvalPred = n.asarray([dval[tuple(i)] \
                for i in izip(n.arange(len(dval),dtype='i4'),indLab[pred])],dtype='f4')
        predNew = (dvalPred>minDval).choose(0,pred)
        #predNew = pred
        return predNew

class Predictor2(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)
        #iDval = n.argsort(dval,-1)
        self.dval.sort(-1)

    def predict(self,minDval):
        dval = self.dval
        pred = self.pred
        predNew = n.logical_and(dval[:,-1] > 0,dval[:,-2] < 0).choose(0,pred)
        return predNew

class Predictor3(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)
        #iDval = n.argsort(dval,-1)
        self.dval.sort(-1)

    def predict(self,minDval):
        dval = self.dval
        pred = self.pred
        predNew = (dval[:,-1] - dval[:,-2] > minDval).choose(0,pred)
        return predNew

class Predictor4(Predictor):

    def __init__(self,**kw):
        Predictor.__init__(self,**kw)
        #iDval = n.argsort(dval,-1)
        self.dval.sort(-1)

    def predict(self,minDval):
        dval = self.dval
        pred = self.pred
        predNew = (n.abs((dval[:,-1] - dval[:,-2])/(dval[:,-1] + dval[:,-2])) > minDval).choose(0,pred)
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
        colToLev = n.zeros((len(labels),len(self.levelNames)),dtype='i4')
        labelToId = self.labelToId
        levels = self.levels
        taxaTree = self.taxaTree
        indLab = self.indLab
        for iCol in xrange(len(labels)):
            taxid = labelToId[labels[iCol]]
            colToLev[iCol,:] = levels.lineageFixedList(node=taxaTree.getNode(taxid),null=0)
        colToLevLab = n.zeros_like(colToLev)
        # levLab[level ind] = { taxid -> label }
        self.levLab = []
        # labLev[level ind] = { label -> taxid }
        self.labLev = []
        # same but { label -> taxa name }
        self.labLevName = []
        for iLev in xrange(colToLev.shape[1]):
            uniqLev = n.unique(colToLev[:,iLev])
            if uniqLev[0] == 0:
                uniqLev = uniqLev[1:]
            levLab = dict(izip(uniqLev,n.arange(1,len(uniqLev)+1,dtype='i4')))
            levLab[0] = 0
            for iCol in xrange(len(colToLev)):
                colToLevLab[iCol,iLev] = levLab[colToLev[iCol,iLev]]
            self.levLab.append(levLab)
            self.labLev.append(dict([ (item[1],item[0]) for item in levLab.items() ]))
            self.labLevName.append(dict([ (item[1],self.getTaxaName(item[0])) 
                for item in levLab.items() ]))
        sh = colToLevLab.shape
        labToLevLab = n.zeros((sh[0]+1,sh[1]),dtype='i4')
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

def loadPredLibLinear(fileName,nTests):
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
            labels = n.fromstring(line[len('labels'):],dtype='i4',sep=' ')
            pred = n.zeros(nTests,dtype='i4')
            dval = n.zeros((nTests,len(labels)),dtype='f4')
            iPred = 0
            for line in inp:
                val = n.fromstring(line,dtype='f4',sep=' ')
                pred[iPred] = int(val[0])
                dval[iPred] = val[1:]
                iPred += 1
                if iPred % 10000 == 0:
                    print "Loaded %s/%s predictions" % (iPred,len(pred))
            assert iPred == nTests, "Number of tests = %s, number of predictions = %s" % (nTests,iPred)
            indLab = n.zeros(labels.max()+1,dtype='i4')
            indLab[:] = 999999999
            for i in xrange(len(labels)): indLab[labels[i]] = i
            res = (pred,labels,indLab,dval)
        else:
            inp.seek(0)
            pred = n.fromfile(inp,dtype='i4',sep='\n')
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
