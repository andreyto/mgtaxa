"""Scripts for phage/host assignment"""

from MGT.Phage import *
from MGT.Taxa import *
from MGT.Svm import *

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

