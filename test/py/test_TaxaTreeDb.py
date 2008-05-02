from MGT.Taxa import *

from MGT.Config import options

storeDump = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile)

db = createDbSql()

tableSfx = "tmp_test"
tableSfxCmp = "tmp_test_cmp"

storeDb = NodeStorageDb(db=db,tableSfx=tableSfx)

taxaTree = TaxaTreeDb(storage=storeDump,db=db,tableSfx=tableSfx)
taxaTree.loadSeqLen()

#for node in taxaTree.iterDepthTop():
#    print node

taxaLevels = TaxaLevels()
taxaLevels.setLevels(taxaTree)

storeDb.save(taxaTree)

#sys.exit(0)

taxaTreeCmp = TaxaTreeDb(storage=storeDb,db=db,tableSfx=tableSfx)

storeDbCmp = NodeStorageDb(db=db,tableSfx=tableSfxCmp)

storeDbCmp.save(taxaTreeCmp)

#for node in taxaTree.getNodesIter():
#    print taxaLevels.lineageKeys(node), node
