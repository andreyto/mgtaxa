"""Emulate hierachical database-like interface using directiry tree as storage"""

from MGT.Svm import *
from MGT.App import *
from MGT.ClassifierApp import *
from MGT.ParamScanApp import *
from MGT.SeqFeaturesApp import *
import glob

class DirStore:

    metaDirName = ".dirstore"
    dumpName = "store.pkl"
    stemName = "all"
    objExt = ".pkl"

    @classmethod
    def open(klass,path,mode="a",**kw):
        dPath = klass._dumpPath(path)
        if mode == "w":
            mode = "c" #currently we just re-create the store state files
        if mode == "a":
            if not os.path.isfile(dPath):
                mode = "c"
            else:
                mode = "r"
        if mode == "r":
            o = loadObj(dPath)
            o.path = path
            # fix path so that we can safely move store trees on disk
        elif mode == "c":
            o = klass(path=path,**kw)
            o.save()
        else:
            raise ValueError("%s - unknown mode value" % mode)
        return o

    def __init__(self,path,opt=None):
        self.path = path
        if opt is None:
            opt = Struct()
        self.opt = opt

    def save(self):
        mDir = self._metaPath(self.path)
        if not os.path.isdir(mDir):
            os.makedirs(mDir)
        dPath = self._dumpPath(self.path)
        dumpObj(self,dPath)

    @classmethod
    def _metaPath(klass,path):
        return pjoin(path,klass.metaDirName)

    @classmethod
    def _dumpPath(klass,path):
        return pjoin(klass._metaPath(path),klass.dumpName)

    @classmethod
    def isStore(klass,path):
        return os.path.isdir(klass._metaPath(path))

    def getFilePath(self,name):
        return pjoin(self.path,name)
    
    def getObjPath(self,name):
        assert name != self.metaDirName
        return self.getFilePath(name+self.objExt)

    def getStorePath(self,name):
        assert name != self.metaDirName
        return self.getFilePath(name)
    
    def getPath(self):
        return self.path

    def getName(self):
        return os.path.basename(self.path)

    def hasFile(self,name):
        return os.exists(self.getFilePath(name))
    
    def hasObj(self,name):
        return os.exists(self.getObjPath(name))

    def hasStore(self,name):
        return self.isStore(self.getStorePath(name))

    def loadObj(self,name):
        objPath = self.getObjPath(name)
        if os.path.isfile(objPath):
            return loadObj(objPath)
        else:
            objPath = self.getStorePath(name)
            if self.isStore(objPath):
                return self.open(path=objPath,mode="r")
            else:
                raise ValueError("Name %s does not point to a valid object" % name)

    def loadStore(self,name):
        return self.open(path=self.getStorePath(name),mode="r")

    def saveObj(self,obj,name):
        if not isinstance(obj,DirStore):
            dumpObj(obj,self.getObjPath(name))
        else:
            obj.path = self.getStorePath(name)
            obj.save()

    def subStore(self,name,klass=None,**kw):
        if klass is None:
            klass = self.__class__
        return klass.open(path=self.getStorePath(name),**kw)

    def objNames(self,pattern="*"):
        for path in glob.iglob(self.getFilePath(pattern+self.objExt)):
            yield os.path.basename(path)

    def storeNames(self,pattern="*"):
        for path in glob.iglob(self.getFilePath(pattern)):
            if self.isStore(path):
                yield os.path.basename(path)

    def objects(self,pattern="*",klass=None):
        for name in self.objNames(pattern=pattern):
            o = self.loadObj(name)
            if klass is None or isinstance(o,klass):
                yield (name,o)

    def stores(self,pattern="*",klass=None):
        if klass is None: #late binding
            klass = DirStore
        for name in self.storeNames(pattern=pattern):
            o = self.loadStore(name)
            if isinstance(o,klass):
                yield o

    def getDefault(self,value,name):
        """Convenience method for getting default values for object names.
        Define a (class) attribute with name 'name'.
        Call this, e.g. o.getDefault(1.,"cutoff") - returns 1., but
        o.getDefault(None,"cutoff") returns value of o.cutoff"""
        if value is None:
            return getattr(self,name)
        else:
            return value

class SampStore(DirStore):

    idMapName = "idmap"
    ## name of default idLabs object
    idLabsName = "idlab"
    ## name of dirstore under which all ParamScanStore's are created
    parScanSupName = "parscan"
    ## name of default ParamScanStore
    parScanName = "1"
    ## name of default sample file
    sampName = "samp"
    

    def getSampFilePath(self,name=None):
        return self.getFilePath(self.getDefault(name,"sampName"))
    
    def fromPreProc(self,inpSamp,preProc):
        data = loadSeqs(inpFile=inpSamp,preProc=preProc)
        saveSeqs(data,self.getSampFilePath())
        self.saveObj(preProc.getIdMap(),self.idMapName)

    def loadIdMap(self):
        return self.loadObj(self.idMapName)

    def loadIds(self):
        return loadSeqsIdDef(self.getSampFilePath())

    def featStore(self,name,**kw):
        return self.subStore(name=name,klass=FeatStore,**kw)

    def featStores(self,pattern="*"):
        for x in self.stores(pattern=pattern,klass=FeatStore):
            yield x

    def makeIdLabs(self,idLabsSrc,name=None,action="idmap",actionOpt={}):
        """Build and save IdLabels using IdLabels object for source samples and a given action.
        @param idLabsSrc IdLabels instance for ids of the source samples (e.g. those passed to byPreProc()).
        @param name if None, default name will be used for saving object, otherwise this name
        @param action a choice of "idmap","idfilt" or "balance"
        @param actionOpt keyword parameters if needed by the action
        @return new IdLabels object (it will also be saved as name)
        After this call, you can use loadIdLabs() to get the IdLabels object again.
        """
        print "Input idLabs counts for action '%s' in store '%s': %s" % (action,self.getName(),idLabsSrc.count(),)
        name = self.getIdLabsName(name)
        if action == "idmap":
            idMap = self.loadIdMap()
            idLabs = idLabsSrc.remapIds(idMap=idMap)
        elif action == "idfilt":
            ids = self.loadIds()
            idLabs = idLabsSrc.selById(ids=ids)
        elif action == "balance":
            idLabs = idLabsSrc.balance(**actionOpt)
        self.saveObj(idLabs,name)
        print "Output idLabs counts for action '%s' in store '%s': %s" % (action,self.getName(),idLabs.count(),)
        return idLabs

    def loadIdLabs(self,name=None):
        name = self.getIdLabsName(name)
        return self.loadObj(name)

    def saveIdLabs(self,idLabs,name=None):
        name = self.getIdLabsName(name)
        return self.saveObj(idLabs,name)

    def getIdLabsPath(self,name=None):
        return self.getObjPath(self.getIdLabsName(name))

    def getIdLabsName(self,name=None):
        if name is None:
            return self.idLabsName
        else:
            return name

    def parScanSupStore(self):
        return self.subStore(name=self.parScanSupName,klass=DirStore)

    def parScan(self,opt,name=None,idLabsName=None):
        store = self.parScanSupStore().subStore(name=self.getDefault(name,"parScanName"),
                klass=ParamScanStore,
                mode="w",
                opt=opt,
                sampStore=self,
                idLabsName=idLabsName)

        return store.run()
    
    def parScanStore(self,name=None):
        return self.parScanSupStore().loadStore(name)



class FeatStore(SampStore):

    def fromOther(self,sampStore):
        opt = self.opt
        ## predefined opt.inSeq should be relative name or None
        opt.inSeq = sampStore.getSampFilePath(opt.get("inSeq",None))
        ## predefined opt.outFeat should be relative name or None
        opt.outFeat = self.getSampFilePath(opt.get("outFeat",None))
        opt.runMode = "inproc"
        opt.cwd = self.getPath()
        app = SeqFeaturesApp(opt=opt)
        app.run()

class ParamScanStore(DirStore):
    
    def __init__(self,path,opt,sampStore=None,idLabsName=None):
        opt = deepcopy(opt)
        if sampStore is not None:
            o = opt.clOpt
            o.inFeat = [ sampStore.getSampFilePath() ]
            if o.get("labels",None) is None:
                o.labels = [sampStore.getIdLabsPath(idLabsName)]
        opt.cwd = path
        DirStore.__init__(self,path,opt=opt)

    def run(self):
        opt = self.opt
        opt.cwd = self.path # in case we moved it after __init__()
        opt.mode = "scatter"
        app = ParamScanApp(opt=opt)
        self.opt = app.getOpt()
        self.save()
        return app.run()

