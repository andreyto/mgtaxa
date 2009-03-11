from MGT.Taxa import *
from MGT.SeqIO import pullNCBISeqIds
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


class SynecApp(App):
    """App-derived class for synecoccus classification"""

    batchDepModes = ("parscan",)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        clOpt = opt.setdefault("clOpt",Struct())
        psOpt = opt.setdefault("psOpt",Struct())
        ftOpt = opt.setdefault("ftOpt",Struct())
        #clOpt.thresh = n.arange(1.,1.01,0.1)
        #clOpt.saveConfMatr = True
        #clOpt.exportConfMatr = True
        clOpt.method = "svm" #"svm" "knn"
        clOpt.kernel = "lin" #"rbf"
        clOpt.knnK = 10
        ClassifierApp.fillWithDefaultOptions(clOpt)
        pgen = ParamGridGen()
        #psOpt.params = pgen.add("C",[0.00001]).grid() #(-10,10,2)(0,1,2)
        #psOpt.params = pgen.add("C",pgen.p2(6,16,2)).add("rbfWidth",pgen.p2(-4,10,2)).grid()
        psOpt.params = pgen.add("C",pgen.p2(-8,16,1)).grid() #(-10,10,2)(0,1,2)
        #psOpt.params = pgen.add("knnK",pgen.lin(1,20,2)).grid()
        #psOpt.params = pgen.add("knnMaxDist",pgen.lin(0.01,0.1,0.01)).grid()
        psOpt.clOpt = clOpt
        #ftOpt.inSeq = "samp.rndGenGcFixed"
        ftOpt.featType = "kmer" #"wdh" "kmerlad" #"kmer"
        if ftOpt.featType == "wdh":
            ftOpt.revCompl = "merge" #"addcol" #"forward"
        elif ftOpt.featType == "kmerlad":
            ftOpt.kmerLen = 9
            ftOpt.revCompl = "addcol" #"forward"
            ftOpt.norm = NORM_POLICY.EXPECT | NORM_POLICY.EU_ROW
        elif ftOpt.featType == "kmer":
            ftOpt.kmerLen = 8
            ftOpt.revCompl = "merge" #"forward"
            ftOpt.norm = NORM_POLICY.FREQ | NORM_POLICY.EU_ROW
        ftOpt.balance = -2 #do not even shuffle (somehow it takes a while) - we do it anyway when making idlabs
    
    def init(self):
        self.orgs = n.asarray(
                [
                    (64471, True,1),
                    (110662, False,1),
                    (84588, False,2),
                    (316279, True,2),
                    (313625, False,1),
                    (32051, False,2),
                    (59931, False,1),
                    (221360, False,2),
                    (221359, False,1),
                    (316278, False,2)
                ],
                dtype=[("taxid","i4"),("coast",bool),("split","i4")]
                )

        self.orgsBg = n.asarray(
                [   (314261,1),
                    (335992,2),
                    (398580,1),
                    (290400,2),
                    (59922,1),
                    (74547,2),
                    (375451,1),
                    (246200,1),
                    (292414,2)
                ],
                dtype=[("taxid","i4"),("split","i4")]
                )
        self.data = Struct() #that will be saveState() and loadState() content
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        self.seqFile = self.store.getSampFilePath()
    
    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "seq":
            inFastaHdr = [ pjoin(options.refSeqDataDir,"microbial.genomic.fna.hdr") ]
            inFasta    = [ pjoin(options.refSeqDataDir,"microbial.genomic.fna") ]
            self.collectIds(inFastaFiles=inFastaHdr)
            self.checkSeqRecs()
            self.saveState("seqid")
            self.pullSeq(inFastaFiles=inFasta)
        elif opt.mode == "shred":
            self.loadState("seqid")
            for sampLen in opt.sampLens:
                self.prepShreds(sampLen=sampLen)
        elif opt.mode == "feat":
            self.loadState("seqid")
            for sampLen in opt.sampLens:
                self.prepFeat(sampLen=sampLen)
        elif opt.mode == "idlabs":
            self.loadState("seqid")
            self.makeIdLabs()
        elif opt.mode == "parscan":
            self.loadState("seqid")
            return self.parScan(parScanName=opt.get("parScanName",1))

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTreeNew()
        return self.taxaTree

    def initDbMem(self):
        taxa = []
        for rec in self.orgs:
            t = Struct(taxid=rec['taxid'],coast=rec['coast'],split=rec['split'],bg=False)
            taxa.append(t)
        for rec in self.orgsBg:
            t = Struct(taxid=rec['taxid'],split=rec['split'],bg=True)
            taxa.append(t)
        self.data.taxa = taxa

    def collectIds(self,inFastaFiles):
        giToTaxa = loadGiTaxBinNew()
        taxaTree = self.getTaxaTree()
        taxMis = mapFastaRecordsToTaxaTree(inSeqs=inFastaFiles,taxaTree=taxaTree,giToTaxa=giToTaxa,
                storeHeader=True,storeSeqLen=True)
        self.initDbMem()
        taxa = self.data.taxa
        seqRecs = []
        for t in taxa:
            taxid = t.taxid
            node = taxaTree.getNode(taxid)
            for nodeSeq in node.seq:
                if not "plasmid" in nodeSeq.header.lower():
                    seqRec = Struct(gi=nodeSeq.gi,id=nodeSeq.gi,taxid=node.id,
                            maxLen=-1,header=nodeSeq.header,seqLen=nodeSeq.get("seqLen",0))
                    seqRecs.append(seqRec)
        self.data.seqRecs = seqRecs

    def mapTaxaToSeqRecs(self):
        m = defdict(list)
        for seqRec in self.data.seqRecs:
            m[seqRec.taxid].append(seqRec)
        return m

    def checkSeqRecs(self):
        t2s = self.mapTaxaToSeqRecs()
        taxaTree = self.getTaxaTree()
        missing = []
        for t in self.data.taxa:
            node = taxaTree.getNode(t.taxid)
            try:
                seqs = t2s[t.taxid]
                print "Taxa %s \nNode %s\nSequences\n%s\n" % (t,node.name,seqs)
            except KeyError:
                print "Taxid %s \nNode %s\nSequences not found\n" % (t,node.name)
                missing.append((t,node.name))
            print "***********************\n"
        if len(missing) > 0:
            print "Taxids w/o sequence: \n%s\n" % (missing,)

    def getStateObjName(self,name):
        return name+'.st'

    def getSampName(self,sampLen):
        return "samp_%i" % sampLen

    def saveState(self,name):
        self.store.saveObj(self.data,self.getStateObjName(name))

    def loadState(self,name):
        self.data = self.store.loadObj(self.getStateObjName(name))

    def pullSeq(self,inFastaFiles):
        pullNCBISeqIds(inSeqs=inFastaFiles,seqIds=self.data.seqRecs,outSeq=self.seqFile)

    def getSampStore(self,sampLen,**kw):
        return self.store.subStore(name=self.getSampName(sampLen),klass=SampStore,**kw)

    def getFeatStore(self,sampLen,featName):
        return self.getSampStore(sampLen).loadStore(name=featName)

    def prepShreds(self,sampLen):
        shredder = LoadSeqPreprocShred(sampLen,sampNum=0,sampOffset=0,makeUniqueId=True)
        sampStore = self.getSampStore(sampLen=sampLen,mode="w")
        sampStore.fromPreProc(inpSamp=self.seqFile,preProc=shredder)

    def prepFeat(self,sampLen):
        sampStore = self.getSampStore(sampLen=sampLen)
        ftOpt = self.opt.ftOpt
        featStore = sampStore.featStore(name=self.opt.featName,opt=ftOpt,mode="w")
        featStore.fromOther(sampStore=sampStore)

    def makeSeqRecIdLabs(self):
        """Build IdLabels object from seqRec descriptions"""
        useTaxids1 = set([
        64471,
        110662,
        #84588,
        316279,
        32051,
        59931
        ])
        useTaxids2 = n.asarray(
        [
        (64471,   'O',   1),
        (110662,  'S',   1),
        (84588,   'S',   2),
        (316279,  'S',   1),
        (313625,  'S',   2),
        (32051,   'O',   2),
        (59931,   'O',   1),
        (221360,  'O',   2), #Agaba
        (221359,  'O',   1), #Agaba
        (69042,   'O',   2) #LA Sound
        ],
        dtype=[("taxid","i4"),("label","S1"),("split","i4")]
        )
        useTaxids3 = n.asarray(
        [
        (64471,   'O',   1),# 
        (110662,  'S',   2),#
        (84588,   'S',   3),
        (316279,  'S',   4),
        (313625,  'S',   5),
        (32051,   'O',   6),
        (59931,   'O',   7),
        (221360,  'O',   8), #Agaba
        (221359,  'O',   9), #Agaba
        (69042,   'O',   10) #LA Sound
        ],
        dtype=[("taxid","i4"),("label","S1"),("split","i4")]
        )
        useTaxids = n.asarray(
        [
        (64471,   'O',   1),# 
        (110662,  'S',   1),#
        (84588,   'S',   2),
        (316279,  'S',   2),
        #(313625,  'S',   2),
        (32051,   'O',   1),
        (59931,   'O',   2),
        #(221360,  'O',   3), #Agaba
        #(221359,  'O',   3), #Agaba
        #(69042,   'O',   3) #LA Sound
        ],
        dtype=[("taxid","i4"),("label","S1"),("split","i4")]
        )
        taxidLab = groupRecArray(useTaxids,keyField="taxid")
        seqRecs = self.data.seqRecs
        taxaToSeq = self.mapTaxaToSeqRecs()
        taxaRecs = self.data.taxa
        idLabs = IdLabels.makeRecords(nrec=len(seqRecs))
        labToName = {1:"bg",2:"o",3:"s"}
        iLab = 0
        for trec in taxaRecs:
            taxid = trec.taxid
            seqs = taxaToSeq[taxid]
            split = trec.split
            if trec.bg:
                lab = 1
                #TMP:
                #split = 0
            else:
                if taxid in taxidLab:
                    txLab = taxidLab[taxid][0]
                    #lab = 2 if txLab["label"]=='O' else 3
                    lab = 2 if trec.coast else 3
                    split = txLab["split"]
                else:
                    lab = 0
            for seq in seqs:
                idLab = idLabs[iLab]
                idLab["id"] = str(seq["id"])
                idLab["label"] = lab
                idLab["split"] = split
                iLab += 1
        assert iLab == len(idLabs)
        idLabs = idLabs[idLabs["label"] > 0]
        idLabs = IdLabels(records=idLabs)
        idLabs.setLabNames(labToName)
        return idLabs

    def makeIdLabs(self,idLabsName=None):
        opt = self.opt
        idLabsSeq = self.makeSeqRecIdLabs()
        print "idLabs counts for seqs: %s" % (idLabsSeq.count(),)
        for sampLen in opt.sampLens:
            sampStore = self.getSampStore(sampLen)
            nameIL = sampStore.getIdLabsName(idLabsName)
            idLabsSamp = sampStore.makeIdLabs(idLabsSrc=idLabsSeq,name=nameIL)
            for featStore in sampStore.featStores():
                if opt.featName == featStore.getName():
                    idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsSamp,name=nameIL+".unbal",action="idfilt")
                    idLabsFeat = featStore.makeIdLabs(idLabsSrc=idLabsFeat,name=nameIL,action="balance")

    def parScan(self,parScanName=None,idLabsName=None):
        opt = self.opt
        psOpt = opt.psOpt
        psOpt.runMode = opt.runMode
        jobs = []
        for sampLen in opt.sampLens:
            sampStore = self.getSampStore(sampLen)
            for featStore in sampStore.featStores():
                if opt.featName == featStore.getName():
                    jobs.extend(featStore.parScan(opt=psOpt,name=parScanName,idLabsName=idLabsName))
        return jobs


def run_Synec():
    opt = Struct()
    opt.sampLens = [5000]
    opt.parScanName = "2"
    opt.featName = "kmer_8" #"kmer_8a" #"kmerlad_4f_nnorm_rndGenGcFixed" #"kmerlad_4f" #"wd_2_10"
    opt.runMode = "batchDep" #"inproc" #"batchDep"
    modes = ["feat","idlabs"]
    #modes = ["idlabs"]
    #modes = ["parscan"]
    #modes = ["feat","idlabs","parscan"] #["shred","idlabs"]
    jobs = []
    for mode in modes:
        opt.mode = mode #"parscan" #"idlabs" #"shred" #"seq"
        app = SynecApp(opt=opt)
        jobs = app.run(depend=jobs)


