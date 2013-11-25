from MGT.ImmApp import *
from MGT.Taxa import *

immStore = ImmStore.open(sys.argv[1])

taxaTree = loadTaxaTree()
for idMod in immStore.iterIds():
    taxid = int(idMod)
    modMeta = dict(
        id = idMod,
        taxid = taxid,
        name = taxaTree.getNode(taxid).name,
        is_leaf = 1,
        seq_db_ids=[idMod]
        )
    immStore.saveMetaDataById(idMod,meta=modMeta)
    print idMod, modMeta["name"]

