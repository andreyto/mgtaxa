#!/usr/bin/env python
from MGT.Common import *
from MGT.SeqDbFasta import *
from MGT.ImmApp import *

seqDbPath = sys.argv[1]
immDbPath = sys.argv[2]

dbSeq = SeqDbFasta.open(path=seqDbPath,mode="r")
taxids = set(dbSeq.getTaxaList())
taxids.add(648174) #SAR86

immStore = ImmStore.open(path=immDbPath)

for idImm in immStore.listIds():
    idImm = int(idImm)
    meta = dict(taxid=idImm,is_leaf=(idImm in taxids))
    immStore.saveMetaDataById(idImm,meta)



