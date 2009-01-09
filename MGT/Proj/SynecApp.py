from MGT.Taxa import *
from MGT.SeqIO import pullNCBISeqIds
from MGT.Svm import *


class DirStore:

    def __init__(self,workDir):
        self.stemName = "all"
        self.workDir = workDir
        makedir(self.workDir)

    def getFileName(self,name):
        return pjoin(self.workDir,name)

    def getWorkDir(self):
        return self.workDir

    def makeSubStore(self,name,fact=None,**kw):
        if fact is None:
            fact = self.__class__
        return fact(workDir=self.getFileName(name),**kw)

class SampStore(DirStore):

    def getSampFileName(self):
        return self.getFileName(self.stemName+".samp")

    def getIdMapFileName(self):
        return self.getFileName(self.stemName+".samp.idmap")

    def fromOther(self,inpSamp,preProc):
        data = loadSeqs(inpFile=inpSamp,preProc=preProc)
        saveSeqs(data,self.getSampFileName())
        dumpObj(preProc.getIdMap(),self.getIdMapFileName())

    def getFeatStore(self,name):
        return self.makeSubStore(name="feat_%s" % name,fact=FeatStore)


class FeatStore(DirStore):

    def getFeatFileName(self):
        return self.getFileName(self.stemName+".feat")

    def fromSamp(self,sampStore):
        cmd = "%s -i %s -o %s" % (options.getPyCmd("seqToFeat"),sampStore.getSampFileName(),self.getFeatFileName())
        run(cmd,cwd=self.getWorkDir())


class Synec(App):
    """App-derived class for synecoccus classification"""

    #batchDepModes = ("scatter",)

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
        self.workDir = "."
        self.nameStem = "all"
        self.seqFile = self.getFileName(self.nameStem+".seq")
        self.data = Struct() #that will be save() and load() content
        self.taxaTree = None #will be lazy-loaded
    
    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "seq":
            inFastaHdr = [ pjoin(options.refSeqDataDir,"microbial.genomic.fna.hdr") ]
            inFasta    = [ pjoin(options.refSeqDataDir,"microbial.genomic.fna") ]
            self.collectIds(inFastaFiles=inFastaHdr)
            self.checkSeqRecs()
            self.save("seqid")
            self.pullSeq(inFastaFiles=inFasta)
        elif opt.mode == "shred":
            self.load("seqid")
            sampLen = opt.sampLen
            self.prepShreds(sampLen=sampLen)
            self.prepFeat(sampLen=sampLen)
        elif opt.mode == "cv":
            

    def getFileName(self,name):
        return pjoin(self.workDir,name)

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

    def getStateFileName(self,name):
        return self.getFileName(self.nameStem+'.'+name+'.st.pkl')

    def getSampDir(self,sampLen):
        return self.getFileName("samp_%i" % sampLen)

    def save(self,name):
        fileName = self.getStateFileName(name)
        dumpObj(self.data,fileName)

    def load(self,name):
        fileName = self.getStateFileName(name)
        self.data = loadObj(fileName)

    def pullSeq(self,inFastaFiles):
        pullNCBISeqIds(inSeqs=inFastaFiles,seqIds=self.data.seqRecs,outSeq=self.seqFile)

    def getSampStore(self,sampLen):
        return SampStore(workDir=self.getSampDir(sampLen))

    def prepShreds(self,sampLen):
        shredder = LoadSeqPreprocShred(sampLen,sampNum=0,sampOffset=0,makeUniqueId=True)
        sampStore = self.getSampStore(sampLen=sampLen)
        sampStore.fromOther(inpSamp=self.seqFile,preProc=shredder)

    def prepFeat(self,sampLen):
        sampStore = self.getSampStore(sampLen=sampLen)
        featStore = sampStore.getFeatStore("wd_2_10")
        featStore.fromSamp(sampStore=sampStore)


def run_Synec():
    app = Synec()
    stage = "shred"
    if stage == "seq":
        inFastaHdr = [ pjoin(options.refSeqDataDir,"microbial.genomic.fna.hdr") ]
        inFasta    = [ pjoin(options.refSeqDataDir,"microbial.genomic.fna") ]
        app.collectIds(inFastaFiles=inFastaHdr)
        app.checkSeqRecs()
        app.save("seqid")
        app.pullSeq(inFastaFiles=inFasta)
    elif stage == "shred":
        app.load("seqid")
        sampLen = 5000
        #app.prepShreds(sampLen=sampLen)
        app.prepFeat(sampLen=sampLen)
    return app

