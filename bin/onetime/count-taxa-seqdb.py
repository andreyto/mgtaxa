#!/usr/bin/env python
from MGT.Common import *
from MGT.SeqDbFasta import *
from MGT.TaxaTreeUtils import *

"""Count taxonomic groups represented in SeqDbFasta object.
Produces an SQLite DB file containing a table with a "Linnean" 
(in a broad sense) lineage for each SeqDb record, as well as 
the aggregate count tables for each taxonomic level."""

refDataDir = "../../atovtchi/shannon_viral_IO_paper/shannon_viral_paper_2011.v.1"
workDir = refDataDir
#dbPath = pjoin(workDir,"refseq-taxa")
dbPath = "seq-db.mic"
#taxonomyDir = pjoin(workDir,"taxonomy")

topTaxids =  micTaxids

dbSeq = SeqDbFasta(path=dbPath)
taxids = dbSeq.getTaxaList()

#taxaTree = loadTaxaTree(ncbiDumpFile=pjoin(taxonomyDir,"nodes.dmp"),
#                ncbiNamesDumpFile=pjoin(taxonomyDir,"names.dmp"))
taxaTree = loadTaxaTree()

topNodes = [ taxaTree.getNode(topTaxid) for topTaxid in topTaxids ]

linwr = LinnWriter(taxaTree=taxaTree)
wr = linwr.newWriter()

for taxid in taxids:
    node = taxaTree.getNode(taxid)
    if sum(( node.isUnder(topNode) for topNode in topNodes )):
        wr.send(dict(taxid=taxid,weight=dbSeq.seqLengths(taxid)["len"].sum()))
        #wr.send(dict(taxid=taxid))
        #rd = dbSeq.fastaReader(taxid)
        #for rec in rd.records():
        #    print "%s\t%s\t%s" % (taxid,rec.seqLen(),rec.header())
        #rd.close()
wr.close()
