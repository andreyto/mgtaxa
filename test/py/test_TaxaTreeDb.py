### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Taxa import *

from MGT.Config import options

storeDump = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile)

#db = createDbSql()
db = DbSqlLite(dbpath="tmp.sqlite",strategy="exclusive_unsafe")

tableSfx = "tmp_test"
tableSfxCmp = "tmp_test_cmp"

storeDb = NodeStorageDb(db=db,tableSfx=tableSfx)

taxaTree = TaxaTree(storage=storeDump)

taxaLevels = TaxaLevels(taxaTree=taxaTree)

storeDbLev = LevelsStorageDb(db=db,tableSfx=tableSfx)

storeDbLev.save(taxaLevels)

storeDb.save(taxaTree)

storeDb.saveName(taxaTree)

for iNode,node in enumerate(taxaTree.iterDepthTop(),1):
    node.tmp_attr = iNode
taxaTree.setTotal(srcAttr="tmp_attr",dstAttr="tmp_attr_tot")

tableAttr = 'tmp_attr_'+tableSfx

storeDb.saveAttributes(tree=taxaTree,
        nameToSql={'tmp_attr':{'type':'bigint','createIndex':True,'name':'tmp_attr1'}},
        table=tableAttr)

storeDb.loadAttribute(tree=taxaTree,
                   name="tmp_attr2",
                   sql="select id,tmp_attr1 from %s" % (tableAttr,),
                   default=0L,
                   setDefault=True,
                   ignoreKeyError=True,
                   typeCast=long)
taxaTree.setTotal(srcAttr="tmp_attr2",dstAttr="tmp_attr2_tot")


storeDb.saveAttributes(tree=taxaTree,
        nameToSql={'tmp_attr2':{'type':'bigint','createIndex':True},
                   'tmp_attr2_tot':{'type':'bigint','name':'tmp_attr2_tot'}},
        table=tableAttr)

attrTotal = db.selectScalar("select sum(tmp_attr2) from %s" % (tableAttr,))
rootTotal = taxaTree.getRootNode().tmp_attr_tot
assert attrTotal == rootTotal, "Mismatch between computed and stored total sequence length: attrTotal = %s, rootTotal = %s" % (attrTotal,rootTotal)

#for node in taxaTree.iterDepthTop():
#    print node


#sys.exit(0)

taxaTreeCmp = TaxaTree(storage=storeDb)

storeDbCmp = NodeStorageDb(db=db,tableSfx=tableSfxCmp)

storeDbCmp.save(taxaTreeCmp)

#for node in taxaTree.getNodesIter():
#    print taxaLevels.lineageKeys(node), node
