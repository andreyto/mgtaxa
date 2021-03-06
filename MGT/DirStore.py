### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Emulate hierachical database-like interface using directiry tree as storage"""

from MGT.Svm import *
from MGT.App import *
from MGT.ParamScanApp import *
from MGT.SeqFeaturesApp import *
import glob, shutil

class DirStore:

    metaDirName = ".dirstore"
    dumpName = "store.pkl"
    stemName = "all"
    objExt = ".pkl"
    metaExt = ".meta.pkl"

    @classmethod
    def open(klass,path,mode="a",**kw):
        """Open existing or create a new store (optionally deleting existing one with the same path).
        @param path filesystem path of the store
        @param mode string with one of 
        @li "w" - will open exisiting store or create new and call save() method
        @li "c" - will erase existing store (if exists) then create the new store and call save()
        @li "r" - will open existing store, will never create a new one, will not call save()
        @li "r+" - will open existing store for updates; store must exist
        @li "a" - if store exists, acts same as "r+", else same as "w"."""
        dPath = klass._dumpPath(path)
        if mode == "c":
            rmdir(path)
            mode = "w" 
        if mode == "a":
            if not os.path.isfile(dPath):
                mode = "w"
            else:
                mode = "r+"
        if mode[0] == "r":
            o = loadObj(dPath)
            # fix path so that we can safely move store trees on disk
            o.path = path
            if mode.endswith("+"):
                hashers_mode = "w"
            else:
                hashers_mode = "r"
            o.path_hasher, o.meta_path_hasher = klass._getPathHashers(path,mode=hashers_mode)
        elif mode == "w":
            kw_c = kw.copy()
            kw_c["save"] = True
            o = klass(path=path,**kw_c)
        else:
            raise ValueError("%s - unknown mode value" % mode)
        return o

    def __init__(self,path,opt=None,save=False):
        """Constructor.
        @param save If False [default], this only creates in-memory instance but does not touch
        the file system in any way. This is an intented behavior, so that
        the user code could create a hierarchy of these objects w/o commiting
        anything to disk (if it only needs to get to the leaves, for instance).
        To commit the object to disk, pass save as True. Alternatively, you ca call save() 
        on this instance, or create the instance with a factory class method open().
        """
        self.path = path
        if opt is None:
            opt = Struct()
        self.opt = opt
        if save:
            self.save()
        else:
            self.path_hasher, self.meta_path_hasher = \
                    self._getPathHashers(self.path,mode="r")

    def __getstate__(self):
        """Pickling support - exclude everything but self.opt and self.path to make saved objects lightweight."""
        return dict(opt=self.opt,path=self.path)

    def save(self):
        mDir = self._metaPath(self.path)
        if not os.path.isdir(mDir):
            os.makedirs(mDir)
        dPath = self._dumpPath(self.path)
        dumpObj(self,dPath)
        self.path_hasher, self.meta_path_hasher = self._getPathHashers(self.path,mode="w")

    def close(self):
        pass

    def setAttr(self,key,val):
        self.opt["key"] = val

    def getAttr(self,key,default):
        return self.opt.get(key,default)

    @classmethod
    def _metaPath(klass,path):
        return pjoin(path,klass.metaDirName)

    @classmethod
    def _dumpPath(klass,path):
        return pjoin(klass._metaPath(path),klass.dumpName)
    
    @classmethod
    def _getPathHashers(klass,path,mode):
        return (PathHasher(path,mode=mode),
                PathHasher(klass._metaPath(path),mode=mode))

    @classmethod
    def isStore(klass,path):
        return os.path.isdir(klass._metaPath(path))

    def getFileMetaPath(self,name):
        return self.meta_path_hasher(name+self.metaExt)

    def delFileMetaData(self,name,ignore=True):
        p = self.getFileMetaPath(name)
        try:
            os.remove(p)
        except OSError:
            if not ignore:
                raise
    
    def hasFileMetaData(self,name):
        p = self.getFileMetaPath(name)
        return os.path.exists(p)
    
    def loadFileMetaData(self,name):
        p = self.getFileMetaPath(name)
        return loadObj(p)
    
    def saveFileMetaData(self,obj,name):
        p = self.getFileMetaPath(name)
        dumpObj(obj,p)
    
    def getFilePath(self,name):
        return self.path_hasher(name)
    
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

    def delName(self,name,ignore=True):
        try:
            os.remove(self.getFilePath(name))
        except OSError:
            if not ignore:
                raise
        self.delFileMetaData(name,ignore=ignore)

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
        return self.open(path=self.getStorePath(name),mode="r+")

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
        for path in self.path_hasher.glob(self.getFilePath(pattern+self.objExt)):
            yield stripSfx(os.path.basename(path),self.objExt)

    def getFilePaths(self,pattern="*"):
        for path in self.path_hasher.glob(self.getFilePath(pattern)):
            yield path
    
    def fileNames(self,pattern="*",sfxStrip=None,iterPaths=None):
        """Iterate over names that optionally filtered by the shell pattern.
        @param pattern Shel-style pattern
        @param sfxStrip Filename suffix to strip before returning the remaining
        part (match is done before stripping)
        @param iterPaths If defined, must be an iterable over file paths - it 
        will be used instead of the content of this store. This option can be
        used to build multi-level filters.
        @return iterator of basenames with suffix stripped"""
        if iterPaths is not None:
            import fnmatch
            def _iter_glob():
                for item in iterPaths:
                    if fnmatch.fnmatch(item,pattern):
                        yield item
        else:
            def _iter_glob():
                for item in self.path_hasher.glob(pattern):
                    yield item
        for path in _iter_glob():
            f = os.path.basename(path)
            if sfxStrip is not None:
                f = stripSfx(f,sfxStrip)
            yield f
    
    def storeNames(self,pattern="*"):
        for path in self.path_hasher.glob(pattern):
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

    def openStream(self,name,*l,**kw):
        """Return an open file object stream for name"""
        return openCompressed(self.getFilePath(name),*l,**kw)

class DirKeyStore(DirStore):
    """DirStore that holds one kind of objects and can address them by ID key"""

    ## Suffix for the keyed objects stored in this store. Derived 
    ## classes can override this attribute
    objSfx = ".obj"
    
    def getFileBaseById(self,id,objSfx=None):
        if objSfx is None:
            objSfx = self.objSfx
        return "%s%s" % (id,objSfx)
    
    def getFilePathById(self,id,objSfx=None):
        return self.getFilePath(self.getFileBaseById(id,objSfx=objSfx))
    
    def loadMetaDataById(self,id,default=None):
        if self.hasMetaDataById(id):
            return self.loadFileMetaData(self.getFileBaseById(id))
        elif default is not None:
            return default
        else:
            raise KeyError(id)

    
    def saveMetaDataById(self,id,meta):
        self.saveFileMetaData(meta,self.getFileBaseById(id))
    
    def hasMetaDataById(self,id):
        return self.hasFileMetaData(self.getFileBaseById(id))
    
    def delById(self,id,ignore=True):
        self.delName(self.getFileBaseById(id),ignore=ignore)

    def iterIds(self,iterPaths=None,objSfx=None):
        """Iterate over object IDs either in this store or from the externally provided iterable"""
        if objSfx is None:
            objSfx = self.objSfx
        return self.fileNames(pattern="*"+objSfx,sfxStrip=objSfx,iterPaths=iterPaths)
    
    def iterMetaData(self,iterPaths=None,objSfx=None):
        """Iterate over metadata either in this store or from the externally provided iterable.
        @return iterator over tuples (object id,metadata object)"""
        for id in self.iterIds(iterPaths=iterPaths,objSfx=objSfx):
            yield (id,self.loadMetaDataById(id))
    
    def dictMetaData(self,iterPaths=None,objSfx=None):
        """Load a dict{id->metadata} either from this store or from the externally provided iterable."""
        return dict(self.iterMetaData(iterPaths=iterPaths,objSfx=objSfx))

    def dictMetaDataByIds(self,ids):
        """Return dict{id->metadata} for a sequence of ids"""
        return dict(( (id,self.loadMetaDataById(id)) for id in ids ))
    
    def listIds(self,iterPaths=None,objSfx=None):
        """List object IDs from either this store or from the externally provided iterable"""
        return list(self.iterIds(iterPaths=iterPaths,objSfx=objSfx))
    
    def openStreamById(self,id,*l,**kw):
        return self.openStream(self.getFileBaseById(id),*l,**kw)

    def writeToStreamByIds(self,ids,out):
        for id in ids:
            with closing(self.openStreamById(id,"r")) as inp:
                shutil.copyfileobj(inp,out,1024**2)

class SampStore(DirKeyStore):

    idMapName = "idmap"
    ## name of default idLabs object
    idLabsName = "idlab"
    ## name of default param scan idLabs object - should hold cross-val splits and exclude test samples
    idLabsParScanName = "idlab.parscan"
    ## name of default test idLabs object - should hold one split for union of parscan data and another for test data
    idLabsTestName = "idlab.test"
    ## name of default full training idLabs object - should hold one split for union of all data with known labels
    idLabsTrainName = "idlab.train"
    ## name of dirstore under which all ParamScanStore's are created
    parScanSupName = "parscan"
    ## name of default ParamScanStore
    parScanName = "1"
    ## name of default sample file
    sampName = "samp"
    ## name of dirstore under which all TestStore's are created
    testSupName = "test"
    ## name of default TestStore
    testName = "1"
    ## name of dirstore under which all TrainStore's are created
    trainSupName = "train"
    ## name of default TrainStore
    trainName = "1"
    ## name of dirstore under which all PredStore's are created
    predictSupName = "pred"
    ## name of default PredictStore
    predictName = "1"
    

    def getSampFilePath(self,name=None):
        return self.getFilePath(self.getDefault(name,"sampName"))

    def loadSamp(self,name=None):
        return loadSeqs(inpFile=self.getSampFilePath(name=name))
    
    def saveSamp(self,data,name=None):
        ##@todo make it generic - not hard-wrired sparse features
        return saveSparseSeqs(data,outFile=self.getSampFilePath(name=name))
    
    def fromPreProc(self,inpSamp,preProc):
        data = loadSeqs(inpFile=inpSamp,preProc=preProc)
        saveSeqs(data,self.getSampFilePath())
        if hasattr(preProc,"getIdMap"):
            self.saveObj(preProc.getIdMap(),self.idMapName)

    def loadIdMap(self):
        return self.loadObj(self.idMapName)

    def loadIds(self):
        return loadSeqsIdDef(self.getSampFilePath())

    def featStore(self,name,**kw):
        return self.subStore(name=name,klass=FeatStore,**kw)

    def makeShredStore(self,name,shredOpt={},**kw):
        store = self.subStore(name=name,mode="c",**kw)
        shredOpt["makeUniqueId"] = True
        store.fromPreProc(inpSamp=self.getSampFilePath(),
                preProc=LoadSeqPreprocShred(**shredOpt))
        store.makeIdLabs(idLabsSrc=self.loadIdLabs(),action="idmap")
        return store

    def exportFasta(self,name,lineLen=None):
        data = loadSeqs(self.getSampFilePath())
        saveSeqs(data,outFile=self.getFilePath(name),format="fasta",lineLen=lineLen)

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
    
    def parScan(self,opt,name=None,idLabsName=None,**kw):
        """Create a new substore under parScanSupStore() and run parameter scan where.
        If a substore with a given name already exists, it will be first erased."""
        if idLabsName is None:
            idLabsName = self.getIdLabsName(self.idLabsParScanName)
        store = self.parScanSupStore().subStore(name=self.getDefault(name,"parScanName"),
                klass=ParamScanStore,
                mode="c",
                opt=opt,
                sampStore=self,
                idLabsName=idLabsName)

        return store.run(**kw)
    
    def parScanStore(self,name=None):
        return self.parScanSupStore().loadStore(name)
    
    def testSupStore(self):
        return self.subStore(name=self.testSupName,klass=DirStore)

    
    def test(self,opt,name=None,idLabsName=None,**kw):
        if idLabsName is None:
            idLabsName = self.getIdLabsName(self.idLabsTestName)
        store = self.testSupStore().subStore(name=self.getDefault(name,"testName"),
                klass=TestStore,
                mode="c",
                opt=opt,
                sampStore=self,
                idLabsName=idLabsName)

        return store.run(**kw)
    
    def testStore(self,name=None):
        return self.testSupStore().loadStore(name)

    def trainSupStore(self):
        return self.subStore(name=self.trainSupName,klass=DirStore)
    
    def train(self,opt,name=None,idLabsName=None,**kw):
        if idLabsName is None:
            idLabsName = self.getIdLabsName(self.idLabsTrainName)
        store = self.trainSupStore().subStore(name=self.getDefault(name,"trainName"),
                klass=TrainStore,
                mode="c",
                opt=opt,
                sampStore=self,
                idLabsName=idLabsName)

        return store.run(**kw)
    
    def trainStore(self,name=None):
        return self.trainSupStore().loadStore(name)
    
    def predictSupStore(self):
        return self.subStore(name=self.predictSupName,klass=DirStore)
    
    def predict(self,opt,name=None,idLabsName=None,**kw):
        store = self.predictSupStore().subStore(name=self.getDefault(name,"predictName"),
                klass=PredictStore,
                mode="c",
                opt=opt,
                sampStore=self,
                idLabsName=idLabsName)

        return store.run(**kw)
    
    def predictStore(self,name=None):
        return self.predictSupStore().loadStore(name=self.getDefault(name,"predictName"))

class FeatStore(SampStore):

    def fromOther(self,sampStore,**kw):
        opt = self.opt
        ## predefined opt.inSeq should be relative name or None
        opt.inSeq = sampStore.getSampFilePath(opt.get("inSeq",None))
        ## predefined opt.outFeat should be relative name or None
        opt.outFeat = self.getSampFilePath(opt.get("outFeat",None))
        opt.runMode = "inproc"
        opt.cwd = self.getPath()
        app = SeqFeaturesApp(opt=opt)
        app.run(**kw)

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

    def run(self,**kw):
        opt = self.opt
        opt.cwd = self.path # in case we moved it after __init__()
        opt.mode = "scatter"
        app = ParamScanApp(opt=opt)
        self.opt = app.getOpt()
        self.save()
        return app.run(**kw)

class TestStore(DirStore):
    
    def __init__(self,path,opt,sampStore=None,idLabsName=None):
        opt = deepcopy(opt)
        if sampStore is not None:
            o = opt
            o.inFeat = [ sampStore.getSampFilePath() ]
            if o.get("labels",None) is None:
                o.labels = [sampStore.getIdLabsPath(idLabsName)]
        opt.cwd = path
        DirStore.__init__(self,path,opt=opt)

    def run(self,**kw):
        from MGT.ClassifierApp import ClassifierApp
        opt = self.opt
        opt.cwd = self.path # in case we moved it after __init__()
        opt.predFile = pjoin(opt.cwd,"pred.pkl")
        opt.perfFile = pjoin(opt.cwd,"perf.pkl")
        jobs = []
        opt.mode = "trainScatter"
        app = ClassifierApp(opt=opt)
        opt = app.getOpt() #get the missing options filled with defaults
        # we need to clean up the models directory otherwise we run a danger
        # of picking up old models that are not trained in this session
        rmrf(self.getFilePath(opt.modelRoot))
        jobs = app.run(**kw)
        opt.mode = "test"
        app = ClassifierApp(opt=opt)
        opt = app.getOpt() #get the missing options filled with defaults
        self.save()
        kw = kw.copy()
        kw["depend"] = jobs
        jobs = app.run(**kw)
        return jobs

class TrainStore(DirStore):
    
    def __init__(self,path,opt,sampStore=None,idLabsName=None):
        opt = deepcopy(opt)
        if sampStore is not None:
            o = opt
            o.inFeat = [ sampStore.getSampFilePath() ]
            if o.get("labels",None) is None:
                o.labels = [sampStore.getIdLabsPath(idLabsName)]
        opt.cwd = path
        DirStore.__init__(self,path,opt=opt)

    def run(self,**kw):
        from MGT.ClassifierApp import ClassifierApp
        opt = self.opt
        opt.cwd = self.path # in case we moved it after __init__()
        jobs = []
        opt.mode = "trainScatter"
        app = ClassifierApp(opt=opt)
        opt = app.getOpt() #get the missing options filled with defaults
        # we need to clean up the models directory otherwise we run a danger
        # of picking up old models that are not trained in this session
        rmrf(self.getFilePath(opt.modelRoot))
        return app.run(**kw)


class PredictStore(DirStore):
    
    def __init__(self,path,opt,sampStore=None,idLabsName=None):
        opt = deepcopy(opt)
        if sampStore is not None:
            o = opt
            o.inFeat = [ sampStore.getSampFilePath() ]
            if o.get("labels",None) is None:
                o.labels = [sampStore.getIdLabsPath(idLabsName)]
        opt.cwd = path
        DirStore.__init__(self,path,opt=opt)

    def run(self,**kw):
        from MGT.ClassifierApp import ClassifierApp
        opt = self.opt
        opt.cwd = self.path # in case we moved it after __init__()
        opt.predFile = pjoin(opt.cwd,"pred.pkl")
        jobs = []
        opt.mode = "predict"
        app = ClassifierApp(opt=opt)
        opt = app.getOpt() #get the missing options filled with defaults
        jobs = app.run(**kw)
        self.save()
        return jobs

