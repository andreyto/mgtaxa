from MGT.Taxa import *
from MGT.SeqIO import pullNCBISeqIds
from MGT.Svm import *
from MGT.App import *
from MGT.ClassifierApp import *
from MGT.ParamScanApp import *
from MGT.SeqFeaturesApp import *

from MGT.DirStore import *

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


