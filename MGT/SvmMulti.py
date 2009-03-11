"""Multiclass SVM classifiers."""

from MGT.Shogun import *
from MGT.Svm import *
from shogun.Features import *
from shogun.Classifier import *
from MGT.Factory import *

def softBinDec(y):
    return numpy.sign(y)*(1-numpy.exp(-numpy.abs(y)))

def srm(wf,c,l,y):
    a = 1-l*y
    return numpy.exp(-(wf + c*numpy.select([a>0],[a],default=0).sum())/(c*len(y)))

def srmBinDec(wf,c,l,y):
    rm = srm(wf,c,l,y)
    sy = softBinDec(y)
    return 

class ModelLinIO:
    
    def getData(self,svm):
        d = Struct()
        d.bias = svm.get_bias()
        d.bias_enabled = svm.get_bias_enabled()
        d.w = svm.get_w()
        return d

    def save(self,out,svm):
        d = self.getData(svm)
        dumpObj(d,out)

    def load(self,inp,svm):
        d = loadObj(inp)
        svm.set_bias(d.bias)
        svm.set_bias_enabled(d.bias_enabled)
        svm.set_w(d.w)

class ModelKernIO:
    
    def getData(self,svm):
        d = Struct()
        d.alphas = svm.get_alphas()
        d.bias = svm.get_bias()
        d.bias_enabled = svm.get_bias_enabled()
        d.svs_ind = svm.get_support_vectors()
        return d

    def save(self,out,svm):
        d = self.getData(svm)
        dumpObj(d,out)

    def load(self,inp,svm):
        d = loadObj(inp)
        svm.create_new_model(len(d.svs_ind))
        svm.set_support_vectors(d.svs_ind)
        svm.set_alphas(d.alphas)
        svm.set_bias(d.bias)
        svm.set_bias_enabled(d.bias_enabled)



class SvmShog(FactoryMixin):

    # FactoryMixin provides a classmethod factory() that accepts the same keywords
    # as the derived's class __init__()

    def classify(self,feat):
        pass

    def train(self,feat,labels):
        pass

    def getModIO(self):
        return self.modIO

    def setModIO(self,modIO):
        self.modIO = modIO
    
    def saveMod(self,out):
        self.getModIO().save(out,self.svm)

    def loadMod(self,inp):
        return self.getModIO().load(inp,self.svm)

class SvmShogLin(SvmShog):

    def __init__(self,C=1.,C2=None):
        #Shogun LibLinear is broken, inverts the binary prediction sign for some datasets
        svm = SVMOcas(SVM_OCAS)
        if C2 is None:
            C2 = C
        svm.set_C(C,C2)
        self.svm = svm
        self.setModIO(ModelLinIO())

    def classify(self,feat):
        svm = self.svm
        nSamp = feat.get_num_vectors()
        lab = Labels(n.zeros(nSamp,dtype='f8'))
        svm.set_features(feat)
        w = svm.get_w()
        n_w = len(w)
        if feat.get_num_features() < n_w:
            feat.set_num_features(n_w)
        return svm.classify(lab).get_labels()

    def train(self,feat,labels):
        svm = self.svm
        epsilon=1e-2
        num_threads=1
        svm.set_features(feat)
        svm.set_labels(Labels(labels))
        #svm=SVMOcas(opt.C, feat, Labels(labBin))
        #svm.set_C(opt.C*float((labBin==-1).sum())/(labBin==1).sum(),opt.C)
        svm.set_epsilon(epsilon)
        svm.parallel.set_num_threads(num_threads)
        svm.set_bias_enabled(True)
        svm.train()

class SvmShogKern(SvmShog):

    def __init__(self,kernel,C=1.,C2=None):
        #Shogun LibLinear is broken, inverts the binary prediction sign for some datasets
        svm = SVMLight() #SVMLight(C,kernel)
        if C2 is None:
            C2 = C
        svm.set_C(C,C2)
        svm.set_kernel(kernel)
        self.svm = svm
        self.setModIO(ModelKernIO())

    def classify(self,feat):
        svm = self.svm
        kernel = svm.get_kernel()
        nSamp = feat.get_num_vectors()
        lab = Labels(n.zeros(nSamp,dtype='f8'))
        kernel.init(kernel.get_lhs(),feat)
        return svm.classify(lab).get_labels()

    def train(self,feat,labels):
        svm = self.svm
        kernel = svm.get_kernel()
        epsilon=1e-2
        num_threads=1
        kernel.init(feat,feat)
        svm.set_labels(Labels(labels))
        svm.set_epsilon(epsilon)
        svm.parallel.set_num_threads(num_threads)
        svm.set_bias_enabled(True)
        svm.train()


class SvmModStore:

    def __init__(self,svmFact=None):
        self.svmFact = svmFact

    def getSvmFact(self):
        return self.svmFact

    def setSvmFact(self,svmFact):
        self.svmFact = svmFact

    def newSvm(self):
        return self.svmFact()

    def loadModels(self,labels=None):
        if labels == None:
            labels = self.getLabels()
        for label in labels:
            yield (label,self.loadMod(label))

    def saveMod(self,label,svm):
        pass

    def loadMod(self,label):
        pass

    def getLabels(self):
        pass

    def getMaxLabel(self):
        labels = self.getLabels()
        return labels[-1]

    def fromOther(self,other,labels=None):
        for label,mod in other.loadModels(labels=labels):
            self.saveMod(label,mod)

class SvmModFileStore(SvmModStore):

    maxLabel = 100000
    def __init__(self,modelRoot,**kw):
        SvmModStore.__init__(self,**kw)
        assert '-' not in modelRoot,"model name cannot contain dash symbol"
        self.modelRoot=modelRoot
        self.numLabDig = int(numpy.ceil(numpy.log10(self.maxLabel))+1)

    def name(self,label):
        labFormat = "%0"+str(self.numLabDig)+"i"
        return os.path.join(self.modelRoot,labFormat % label)

    def saveMod(self,label,svm):
        fileName = self.name(label)
        makeFilePath(fileName)
        svm.saveMod(fileName)

    def loadMod(self,label):
        svm = self.newSvm()
        svm.loadMod(self.name(label))
        return svm

    def getLabels(self):
        return sorted([ int(modFile) for modFile in os.listdir(self.modelRoot) ])


class SvmModMemStore(SvmModStore):

    def __init__(self,*l,**kw):
        SvmModStore.__init__(self,*l,**kw)
        self.svms = {}

    def saveMod(self,label,svm):
        self.svms[label] = svm

    def loadMod(self,label):
        return self.svms[label]

    def getLabels(self):
        return sorted(self.svms.keys())


class SvmOneVsAll:

    def __init__(self,maxLabel=None):
        self.maxLabel = maxLabel
        #keep bin svms in memory by default
        self.modStore = SvmModMemStore() 

    def setFeat(self,feat):
        self.feat = feat

    def setLab(self,lab):
        self.lab = lab
        if self.maxLabel is None and self.lab is not None:
            self.maxLabel = max(self.lab)

    def setSvmStore(self,modStore):
        self.modStore = modStore

    def classifyBin(self):
        modIter = self.modStore.loadModels()
        feat = self.feat
        nSamp = feat.get_num_vectors()
        binLab = numpy.zeros((nSamp,self.maxLabel+1),dtype='f8')
        labSeen = numpy.zeros(binLab.shape[1],dtype=bool)
        for (label,svm) in modIter:
            labPred = svm.classify(feat)
            binLab[:,label] = labPred
            labSeen[label] = True
        #assert not labSeen[0] #why did I need this?
        self.binLab = binLab
        self.labSeen = labSeen
        print "SvmOneVsAll.classifyBin() did not see models for labels: "+','.join(["%s" % l for l in n.where(labSeen == False)[0]])

    def trainBin(self,trainLabel):
        feat = self.feat
        lab = self.lab
        labBin = numpy.select([lab == trainLabel], [1.], default=-1.)
        svm = self.modStore.newSvm()
        print "Training label %i ..." % trainLabel
        svm.train(feat=feat,labels=labBin)
        return svm

    def trainMany(self,trainLabels=None):
        if trainLabels is None:
            trainLabels = numpy.unique(self.lab)
        modStore = self.modStore
        for trainLabel in trainLabels:
            svm = self.trainBin(trainLabel)
            modStore.saveMod(trainLabel,svm)

    def computeSrm(self):
        labels = self.lab
        modIter = self.modStore.loadModels()
        self.srm = numpy.zeros(len(self.labSeen),dtype='f8')
        for (label,svm) in modIter:
            ##FIXME for RBF
            w = svm.get_w()
            wf = numpy.inner(w,w)/2.
            c = svm.get_C1()
            assert self.labSeen[label]
            y = self.binLab[:,label]
            l = numpy.select([labels == label],[1.],default=-1.)
            self.srm[label] = srm(wf,c,l,y)

    def getSrm(self):
        return self.srm

    def setSrm(self,srm):
        self.srm = srm

    def classify(self,thresh=-2000,useSrm=False):
        binLab = self.binLab.copy()
        if useSrm:
            for label in numpy.where(self.labSeen)[0]:
                binLab[:,label] = softBinDec(binLab[:,label]) * self.srm[label]
        binLab[:,-self.labSeen] = -1e10
        labMaxBin = binLab.argmax(1)
        rowInd = numpy.arange(len(binLab))
        pred = numpy.select([binLab[rowInd,labMaxBin] >= thresh], [labMaxBin], default=0)
        return pred

def svmOneVsAllOneStep(feat,lab,opt=Struct(thresh=-200,useSrm=False)):
    svmMul = SvmOneVsAll()
    svmMul.setLab(lab[0])
    svmMul.setFeat(feat[0])
    svmMul.trainMany()
    svmMul.setLab(None)
    svmMul.setFeat(feat[1])
    svmMul.classifyBin()
    labPred = svmMul.classify(thresh=opt.thresh,useSrm=opt.useSrm)
    return labPred
