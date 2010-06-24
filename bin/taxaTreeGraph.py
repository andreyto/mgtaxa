### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Taxa import *

from MGT.Config import options

from optparse import OptionParser


def parseCmdLine():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-o", "--out-graph", dest="outGraph",
                    help="output graph file")
    
    (cmdOptions, cmdArgs) = parser.parse_args()
    cmdOptions.outGraph or parser.error("--out-graph is required")
    return (cmdOptions, cmdArgs)

(cmdOptions, cmdArgs) = parseCmdLine()

store = NodeStorageNcbiDump(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile)

db = createDbSql()

taxaTree = TaxaTreeDb(storage=store,db=db)

taxaTree.loadSeqLen()

#taxaTree.deleteNodesIf(lambda n: n.rank == "species")
taxaTree.deleteNodesIf(lambda n: n.id in (28384,410657,12908) or
(n.getParentId() not in (0,1) and (n.getParent().name.lower().startswith("environmental") or n.getParent().name.lower().startswith("unclassified"))))
taxaTree.reindex()

#taxaTree = TaxaTree(storage=store)

#for node in taxaTree.iterDepthTop():
#    print node


taxaLevels = TaxaLevels()
taxaLevels.setLevels(taxaTree)
taxaLevels.reduceNodes(taxaTree)
taxaTree.deleteNodesIf(lambda n: n.getParentId()  == viralRootTaxid and n.rank != "family")
taxaTree.reindex()

#for node in taxaTree.getNodesIter():
#    print node

#storeOut = NodeStorageNewick(fileName="test_TaxaTree.newick",
#    labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))

storeOut = NodeStorageHypView(fileName=cmdOptions.outGraph+".hv3",
    labeler=lambda n: "%s_%s_%s_N_%s_S_%s" % (n.id,n.rank,n.name[:10],n.numChildren(),n.seq_len_tot))

storeOut.save(taxaTree)

storeOut = NodeStorageLibSea(fileName=cmdOptions.outGraph+".libsea",
   labeler=lambda n: "%s_%s_%s" % (n.name,n.rank,n.id))
        
storeOut.save(taxaTree)
