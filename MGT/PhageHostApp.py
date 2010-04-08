### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application for phage/host assignment"""

from MGT.PhageHostDb import *
from MGT.Taxa import *
from MGT.Svm import *
from MGT.App import *
from MGT.DirStore import *
import UUID
from MGT.SeqImportApp import *

class PhageHostApp(App):
    """App-derived class for phage-host assignment"""

    #batchDepModes = ("blastcr","pilercr")

    maxSeqIdLen = UUID.maxIdLen
    maxSeqPartIdLen = UUID.maxIdLen

    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        optChoicesMode = ("build-db-ph","sel-db-ph-pairs","shred-vir-test","predict")
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=optChoicesMode,
            dest="mode",default="predict",help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            make_option(None, "--db-gb-inp",
            action="store", type="string",dest="dbGbInp",help="Input viral GenBank file for phage-host DB construction"),
            make_option(None, "--db-ph",
            action="store", type="string",default="ph",dest="dbPh",help="DB of known phage-host pairs created and used here"),
            make_option(None, "--db-seq",
            action="append", type="string",dest="dbSeqInp",help="Input DB (e.g. RefSeq) FASTA sequence files for both microbes "+\
                    "and viruses referenced by phage-host pairs DB. It can be a subset of superset of p-h DB."),
            make_option(None, "--db-ph-pairs",
            action="store", type="string",default="ph-pairs",dest="dbPhPairs",help="Stratified phage-host pairs from DB phage-host"),
            make_option(None, "--db-ph-pairs-seq-ids",
            action="store", type="string",default="ph-pairs-seq-ids",dest="dbPhPairsSeqIds",help="Sequence IDs for db-ph-pairs - "+\
                    "intermediate file used to pull all sequences from the original FASTA files into db-ph-pairs-seq"),
            make_option(None, "--db-ph-pairs-seq",
            action="store", type="string",default="ph-pairs-seq",dest="dbPhPairsSeq",help="Sequences for db-ph-pairs in our "+\
                    "internal format"),
            make_option(None, "--shred-size-vir-test",
            action="store", type="int",default=400,dest="shredSizeVirTest",help="Shred DB viral samples into this size for testing"),
        ]
        return Struct(usage = "Build phage-host association database, classifiers, and perform training, validation and de novo classification.\n"+\
                "%prog [options]",option_list=option_list)

    
    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        refSeqDir = globals()["options"].refSeqDataDir
        opt.setIfUndef("dbGbInp",pjoin(refSeqDir,"viral.genomic.gbff.gz"))
        opt.setIfUndef("dbSeqInp",[ pjoin(refSeqDir,div+".genomic.fna.gz") for div in ("microbial","viral")])
   

    def init(self):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        self.tmpDir = self.store.getFilePath("tmp")
        makedir(self.tmpDir)
   

    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "build-db-ph":
            return self.buildDbPhageHost()
        elif opt.mode == "sel-db-ph-pairs":
            return self.selectPhageHostPairs()
        elif opt.mode == "shred-vir-test":
            return self.shredVirTest()
        elif opt.mode == "predict":
            return self.predict()
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree


    def loadPhageHostDb(self):
        """Load and map Phage-Host DB into the taxaTree"""
        loadHosts(inpFile=self.store.getObjPath(self.opt.dbPh),taxaTree=self.getTaxaTree())
    
    def buildDbPhageHost(self):
        """Collect and save data about phage-host associations in RefSeq"""
        opt = self.opt
        taxaTree = self.getTaxaTree()
        vhParser = VirHostParser(taxaTree=taxaTree)
        vhParser.assignHostsAndSave(gbFile=opt.dbGbInp,outFile=self.store.getObjPath(opt.dbPh))

    def selectPhageHostPairsTmp(self):
        """Pick pairs of phages and hosts from p-h DB stratified at a genus level"""
        opt = self.opt

        giToTaxa = loadGiTaxBin()
        taxaTree = self.getTaxaTree()


        #taxaTree.getRootNode().setIsUnderUnclass()

        mapFastaRecordsToTaxaTree(inSeqs=opt.dbSeqInp,taxaTree=taxaTree,giToTaxa=giToTaxa)
        self.loadPhageHostDb()

        seqPicker = PhageHostSeqPicker(taxaTree=taxaTree)
        seqPicker.pickSeqHosts()
        seqPicker.groupSeqHosts()
        seqPicker.printGroupSeqHosts()

        seqPicker.pickPairs(maxMicSpe=1,maxMicSeq=1,maxVir=1,maxVirSeq=1)
        seqPicker.checkPairs(giToTaxa=giToTaxa)

        dbPhPairsFile = self.store.getObjPath(opt.dbPhPairs)
        seqPicker.save(dbPhPairsFile)
        seqPicker.load(dbPhPairsFile)
        seqPicker.checkPairs(giToTaxa=giToTaxa)
        seqPicker.saveSeqIds(self.store.getObjPath(opt.dbPhPairsSeqIds))
        
        # import all selected FASTA sequences into internal format

        siOpt = Struct()
        siOpt.runMode = "inproc"
        siOpt.inSeq = opt.dbSeqInp
        siOpt.outFeat = opt.dbPhPairsSeq
        siOpt.inSeqIds = self.store.getObjPath(opt.dbPhPairsSeqIds)
        siOpt.inFormat = "ncbi"
        seqImpApp = SeqImportApp(opt=siOpt)
        seqImpApp.run()
    
    def selectPhageHostPairs(self):
        """Pick pairs of phages and hosts from p-h DB stratified at a genus level"""
        opt = self.opt
        taxaTree = self.getTaxaTree()
        seqPicker = PhageHostSeqPicker(taxaTree=taxaTree)
        dbPhPairsFile = self.store.getObjPath(opt.dbPhPairs)
        seqPicker.load(dbPhPairsFile)
        for div in ("mic","vir","all"):
            self.store.saveIdLabs(seqPicker.makeIdLabsPicks(div=div),name="idlab."+div)

        self.exportVirSamples()

    def getVirSampStoreName(self):
        return "samp.v"
    
    def getVirSampStoreShredsName(self,shredSize):
        return "samp.v.%s" % shredSize

    def getVirSampStore(self):
        return self.store.subStore(self.getVirSampStoreName(),klass=SampStore)

    def getVirSampStoreShreds(self,shredSize):
        return self.getVirSampStore().subStore(self.getVirSampStoreShredsName(shredSize),klass=SampStore)
    
    def exportVirSamples(self):
        opt = self.opt
        sampStore = self.getVirSampStore()
        idLabs = self.store.loadIdLabs("idlab.vir")
        sampStore.fromPreProc(inpSamp=opt.dbPhPairsSeq,
                preProc=LoadSeqPreprocIdFilter(idToLab=idLabs.getIdToLab()))
        sampStore.makeIdLabs(idLabsSrc=idLabs,action="idfilt")
   
    def exportVirSamplesShreds(self,shredSize):
        opt = self.opt
        sampStore = self.getVirSampStore()
        shredStore = sampStore.makeShredStore(name=self.getVirSampStoreShredsName(shredSize),
                shredOpt=dict(sampLen=shredSize,sampNum=10)) #TMP
        shredStore.exportFasta(name="samp.fas",lineLen=1000)

    def shredVirTest(self):
        opt = self.opt
        self.exportVirSamplesShreds(shredSize=opt.shredSizeVirTest)

    def predict(self):
        opt = self.opt
        groupRanks = ("order","family","genus","subgenus","species")
        taxaTree = self.getTaxaTree()
        taxaTree.setMaxSubtreeRank()
        seqPicker = PhageHostSeqPicker(taxaTree=taxaTree)
        dbPhPairsFile = self.store.getObjPath(opt.dbPhPairs)
        seqPicker.load(dbPhPairsFile)
        virToHost = seqPicker.seqVirHostPicks()
        groupRanks = []
        inp = open("/usr/local/projects/GOSII/atovtchi/phymm/results.01.phymm____ph_samp_v_samp_v_400_samp_fas.txt",'r')
        inp.next()
        for line in inp:
            try:
                fields = line.strip().split('\t')
                taxid = int(float(fields[0].split('|')[1]))
                node = taxaTree.getNode(taxid)
                fields[0] = node.name
                predNameGen = fields[3]
                predNode = taxaTree.searchName(predNameGen)
                if len(predNode):
                    predNode = predNode[0]
                    hostNode = virToHost[node][0]
                    groupNode = hostNode.lcsNode(predNode)
                    #groupNode = hostNode.findLeftRankInLineage(groupRanks)
                    if groupNode:
                        groupName = groupNode.name
                        groupRank = groupNode.rank_max
                    else:
                        groupName = "None"
                        groupRank = "no_rank"
                    score = float(fields[2])
                    if score >= - 500:
                        print '\t|\t'.join([groupRank,groupName]+fields)
                        groupRanks.append(groupRank)
            except BaseException, msg:
                print "Exception caught for this line, skipping: %s \t %s" % (msg,line)
        inp.close()
        print binCount(groupRanks,format="list")

class VirHostClassifierScr:

    def __init__(self,virHostsFile=None):
        taxaTree = loadTaxaTreeNew()
        self.taxaTree = taxaTree
        self.micNodes = [ taxaTree.getNode(id) for id in micTaxids ]
        self.virNode = taxaTree.getNode(virTaxid)
        if virHostsFile is not None:
            self.loadVirHosts(virHostsFile)

    def loadVirHosts(self,fileName):
        phost = PhageHostSeqPicker(taxaTree=self.taxaTree)
        phost.load(fileName)
        self.phost = phost
        self.virHosts = phost.seqVirHosts()

    def splitFeatFileIntoVirHost(self,inpFiles,outName):
        data = loadSeqsMany(inpFiles)
        ids = set(data["id"])
        idToNode = {}
        for id in ids:
            idToNode[id] = self.getNCBINode(id)
        virNode = self.virNode
        isVir = n.asarray([ idToNode[id].isSubnode(virNode) for id in data["id"] ],dtype=bool)
        dataVir = data[isVir]
        dataMic = data[isVir!=True]
        labRenum = setLabelsFromIds(dataMic)
        micIdToLab = labRenum.toNew()
        virHostPicks = self.phost.seqVirHostPicks()
        idHostNoSamp = set()
        for rec in dataVir:
            idHost = "%s" % virHostPicks[idToNode[rec["id"]]][0].id
            try:
                rec["label"] = micIdToLab[idHost]
            except KeyError:
                idHostNoSamp.add(idHost)
                rec["label"] = 0
        if len(idHostNoSamp) > 0:
            print "Some hosts have no samples: %s" % (sorted(idHostNoSamp),)
        saveSeqs(dataVir,outName+".v.svm")
        saveSeqs(dataMic,outName+".m.svm")
        dumpObj(labRenum,outName+".labren")

    def getNCBINode(self,idSamp):
        return self.taxaTree.getNode(int(idSamp))



def run_splitFeatFileIntoVirHost():
    vhc = VirHostClassifierScr(virHostsFile="../viralHostsPick.pkl")
    vhc.splitFeatFileIntoVirHost(inpFiles=["all.svm"],outName="all")


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(PhageHostApp)
