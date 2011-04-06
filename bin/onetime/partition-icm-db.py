#!/usr/bin/env python
from MGT.Common import *
from MGT.ImmApp import *

"""Create a bacterial-only Partition of RefSeq ICM database using symlinks.
Creates one directory with symlinks to bac ICMs"""

workDir = os.environ["GOSII_WORK"]
immDbPath = pjoin(workDir,"icm-refseq")
newStorePath = immDbPath+".bac"
newStoreTopTaxid =  bacTaxid

immStore = ImmStore.open(path=immDbPath)
immIds = immStore.listImmIds()

#assume idImm are str(taxids):
immIds = n.asarray(immIds,dtype=int)

taxaTree = loadTaxaTree()

immStoreNew = ImmStore.open(path=immDbPathNew)

topNode = taxaTree.getNode(newStoreTopTaxid)

newFullStorePathRel = "../"+os.path.basename(newFullStorePath)

for immId in immIds:
    node = taxaTree.getNode(immId)
    if node.isUnder(topNode):
        inpFile = immStore.getImmPath(immId)
        outFile = immStoreNew.getImmPath(immId)
        targetPath = pjoin(newFullStorePathRel,os.path.basename(inpFile))
        editSymlink(targetPath,outFile)

