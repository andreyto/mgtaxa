from MGT.Taxa import *

from MGT.Config import options

storeDump = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile)

db = createDbSql()

tableSfx = "tmp_test"
tableSfxCmp = "tmp_test_cmp"

storeDb = NodeStorageDb(db=db,tableSfx=tableSfx)

taxaTree = TaxaTree(storage=storeDump)
storeDb.loadSeqLen(taxaTree)

tableAttr = 'tmp_attr_'+tableSfx

storeDb.saveAttributes(tree=taxaTree,
        nameToSql={'seq_len':{'type':'bigint','createIndex':True},
                   'seq_len_tot':{'type':'bigint','name':'seq_len_tot'}},
        table=tableAttr)

attrTotal = db.selectScalar("select sum(seq_len) from %s" % tableAttr)
rootTotal = taxaTree.getRootNode().seq_len_tot
assert attrTotal == rootTotal, "Mismatch between computed and stored total sequence length: attrTotal = %s, rootTotal = %s" % (attrTotal,rootTotal)

#for node in taxaTree.iterDepthTop():
#    print node

taxaLevels = TaxaLevels()
taxaLevels.setLevels(taxaTree)

storeDb.save(taxaTree)

#sys.exit(0)

taxaTreeCmp = TaxaTree(storage=storeDb)

storeDbCmp = NodeStorageDb(db=db,tableSfx=tableSfxCmp)

storeDbCmp.save(taxaTreeCmp)

#for node in taxaTree.getNodesIter():
#    print taxaLevels.lineageKeys(node), node
