#!/usr/bin/env python
from MGT.Common import *
from MGT.SeqDbFasta import *
from MGT.Taxa import *
from MGT.Functional import *

"""Count taxonomic groups represented in SeqDbFasta object.
Produces an SQLite DB file containing a table with a "Linnean" 
(in a broad sense) lineage for each SeqDb record, as well as 
the aggregate count tables for each taxonomic level."""


class LinnWriter:
    """Writer of a Linnean lineage for a sequence of taxids"""
    
    def __init__(self,taxaTree=None,taxaLevels=None):
        self.taxaTree = taxaTree
        self.taxaLevels = taxaLevels

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree(pklFile=self.opt.taxaTreePkl)
        return self.taxaTree

    def getTaxaLevels(self):
        if self.taxaLevels is None:
            #that assigns "level" and "idlevel" attributes to TaxaTree nodes
            self.taxaLevels = TaxaLevels(self.getTaxaTree())
        return self.taxaLevels

    @coroutine
    def newWriter(self,csvPath="taxa_linn.csv",
            dbPath="taxa_linn.sqlt",
            dbTable="taxa_linn",
            otherFields=list()):
        taxaTree = self.getTaxaTree()
        taxaLevels = self.getTaxaLevels()
        levNames = taxaLevels.getLevelNames("ascend")
        flds = list(otherFields)
        for reqfld in ("weight","rank","name","taxid"):
            if reqfld not in flds:
                flds.insert(0,reqfld)
        for lev in levNames:
            flds += ["taxid_"+lev,"name_"+lev]
        out = open(csvPath,"w")
        w = csv.DictWriter(out, fieldnames=flds, restval='null',dialect='excel-tab')
        w.writerow(dict([(fld,fld) for fld in flds]))
        #w/o handling GeneratorExit exception, the exception
        #is silently handled by the interpreter in the surrounding
        #scope after it is emitted by the 'yield' clause - in
        #other words, the code after the 'while' loop is never
        #called
        try:
            while True:
                rowDict = dict((yield))
                node = taxaTree.getNode(rowDict["taxid"])
                rowDict["name"] = node.name
                rowDict["rank"] = node.linn_level
                rowDict.setdefault("weight",1.0)
                lin = taxaLevels.lineage(node,withUnclass=False)
                for ln in lin:
                    rowDict["name_"+ln.level] = ln.name
                    rowDict["taxid_"+ln.level] = ln.id
                w.writerow(rowDict)
        except GeneratorExit:
            out.close()
            db = DbSqlLite(dbpath=dbPath)
            db.createTableFromCsv(name=dbTable,
                    csvFile=csvPath,
                    hasHeader=True,
                    fieldsMap={"weight":SqlField(type="real")})
            self._sqlReport(db=db,dbTable=dbTable,levNames=levNames)
            db.close()

    def _sqlReport(self,db,dbTable,levNames):
        for levName in levNames:
            fldGrp = "name_"+levName
            db.createTableAs(dbTable+"_grp_"+levName,
                    """\
                    select %(fldGrp)s,
                    sum(weight) as weight
                    from %(dbTable)s
                    group by %(fldGrp)s
                    """ % dict(fldGrp=fldGrp,dbTable=dbTable)
                    )

workDir = pjoin(os.environ["GOSII_WORK"],"shannon_viral_paper_2011")
dbPath = pjoin(workDir,"refseq-taxa")
taxonomyDir = pjoin(workDir,"taxonomy")

topTaxids =  micTaxids

dbSeq = SeqDbFasta(path=dbPath)
taxids = dbSeq.getTaxaList()

taxaTree = loadTaxaTree(ncbiDumpFile=pjoin(taxonomyDir,"nodes.dmp"),
                ncbiNamesDumpFile=pjoin(taxonomyDir,"names.dmp"))

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

