"""Benchmark support for ImmClassifierApp"""

from MGT.DirStore import *
from MGT.SeqDbFasta import *
from MGT.FastaIO import *
from Taxa import *

class ImmClassifierBenchmark(SeqDbFasta):
    """Class that creates shredded benchmarking fragments from SeqDbFasta object.
    The fragment store will use the same SeqDbIds as the source store.
    Provides a method to map fragment IDs to score IDs based on score 
    meta-data"""

    fastaSfx = ".fna"
    
    def __init__(self,*l,**kw):
        DirStore.__init__(self,*l,**kw)
        self.seqDb = None
    
    def getSeqDbIds(self):
        return self.seqDb.dictMetaData()
    
    def getSeqDbIdsImm(self,immDbs):
        immIds = set()
        for immDb in immDbs:
            immIds |= set(immDb.listSeqDbIds())
        return immIds

    def selectIdsDb(self,
            immDbs=None,
            dbBenchTaxidsExcludeTrees=None,
            dbBenchTaxidsInclude=None,
            taxaTree=None):
        """Pick those IDs from SeqDb that will be used to build the benchmark.
        @param immDbs sequence of ImmDb instances - this is used to pick only those entries
        in SeqDb that also have models built against them. If None, no
        filtering by model IDs will be done."""
        idsMeta = self.getSeqDbIds()
        ids = set(idsMeta)
        if immDbs:
            ids &= self.getSeqDbIdsImm(immDbs)
        if dbBenchTaxidsExcludeTrees is not None:
            subTreesExcl = [ taxaTree.getNode(taxid) for taxid \
                    in set(opt.dbBenchTaxidsExcludeTrees) ]
            ids = [ id for id in ids \
                    if not taxaTree.getNode(idsMeta[id]["taxid"]).isUnderAny(subTreesExcl) ] 
        if dbBenchTaxidsInclude is not None:
            with closing(openCompressed(opt.dbBenchTaxidsInclude,'r')) as inp:
                taxidsIncl = set([ int(line.strip()) for line in inp if len(line.strip())>0 ])
            ids = [ id for id in ids if idsMeta[id]["taxid"] in taxidsIncl ]
        return ids

    def makeSample(self,idDb,fragLen,fragCountMax):
        """Make a sample FASTA file for a given SeqDb ID.
        @param idDb Sequences with this SeqDbId will be pulled from
        seqDb attribute, shredded and saved in this object under
        the same SeqDbId"""
        self.shredFasta(idDb=idDb,fragLen=fragLen,
                fragCountMax=fragCountMax)
        self.finById(idDb)

    def catSamples(self,outFile,idsDb=None):
        """Concatenate all sample files for a list of IDs into a single file.
        @note This cats entire files, and thus relies on the newline symbol
        being present at the end of every file. The FastaIO module that is
        used to create individual files always terminates with a newline. """
        import shutil
        
        if idsDb is None:
            idsDb = self.listDbIds()
        if not hasattr(outFile,"write"):
            outFile = openCompressed(outFile,"w")
            outClose = True
        else:
            outClose = False
        for idDb in idsDb:
            inpFile = self.openStreamById(idDb)
            shutil.copyfileobj(inpFile,outFile,1024*1024)
            inpFile.close()
        if outClose:
            outFile.close()

    def shredFasta(self,idDb,fragLen,
            fragCountMax,lineLen=80,outMode="w"):
        """Shred each record in multi-FASTA file into multiple records of fixed size"""
        ##@todo See if we can avoid saving empty files (for very short input records)
        from MGT.FeatIO import LoadSeqPreprocShred
        seqDb = self.seqDb
        seqLenTot = seqDb.seqLengths(idDb)["len"].sum()
        seqLenRatio = (fragCountMax*fragLen)/float(seqLenTot)

        if seqLenRatio < 1.:
            #@todo Pick coords on a virtual concatenation of all sequences, otherwise
            #it will never be quite right.
            sampNum = lambda lab,seq,id: int(rndRound(len(seq)*seqLenRatio/fragLen))
        else:
            sampNum = -1 #all samples
        inpSeq = seqDb.fastaReader(idDb)
        outSeq = self.fastaWriterUncompr(idDb,lineLen=lineLen,mode=outMode)
        shredder = LoadSeqPreprocShred(sampLen=fragLen,
                sampNum=sampNum,
                makeUniqueId=False,
                sortByStarts=True)
        catLen = min(max(seqLenTot/(fragCountMax/10),fragLen*fragCountMax*10),10**8)
        seqCat = []
        seqCatLen = 0
        seqCatStart = 0
        indFrag = 0

        def _out_frags(seqCat,ind,seqCatStart,
                shredder=shredder,idDb=idDb,outSeq=outSeq):
            seqCat = "N".join(seqCat)
            labFr,seqFr,idFr = shredder(0,seqCat,idDb)
            startsFr = shredder.getLastSampStarts()
            for iF,sF in enumerate(seqFr):
                stF = startsFr[iF]
                idF = "%s %s %s pos=%s nincat=%s len=%s" % (genId(),idDb,ind+iF,seqCatStart+stF,iF,len(sF))
                outSeq.record(header=idF,sequence=sF)
            return ind+len(seqFr)
        
        for rec in inpSeq.records():
            hdr = rec.header()
            id = rec.getId()
            seq = rec.sequence()
            seqCat.append(seq)
            seqCatLen += len(seq)
            if seqCatLen >= catLen:
                indFrag = _out_frags(seqCat,indFrag,seqCatStart)
                seqCat = []
                seqCatStart += seqCatLen
                seqCatLen = 0

        if seqCatLen > 0:
            indFrag = _out_frags(seqCat,indFrag,seqCatStart)

        inpSeq.close()
        outSeq.close()
        return indFrag

    def mapIdFragToIdScore(self,idScoreMeta,idSeqDbToIdScoreRemapping=None):
        """Return a mapping from fragment IDs to correct model (score) ID.
        This is used by benchmarking code to get the expected correct prediction
        for each testing fragment.
        @param idScoreMeta mapping of meta data for each model that was used in
        making predictions for the fragments
        @param idSeqDbToIdScoreRemapping optional map to convert each idSeqDb to idScore
        @return dict(idFrag -> idScore)
        """
        if idSeqDbToIdScoreRemapping is None:
            idSeqDbToIdScoreRemapping = dict()
        idSeqDbToIdScore = dict()
        #Map IdSeqDb from score meta-data to score ID
        #We expcect that the self object has benchmarking fragments
        #grouped under the same IdSeqDbs due to the way
        #it was constructed unless idSeqDbToIdScoreRemapping
        #is not None
        for (idScore,scoreMeta) in idScoreMeta.items():
            for idSeqDb in scoreMeta["seq_db_ids"]:
                idSeqDbToIdScore[idSeqDb] = idScore
        
        idFragToIdScore = dict()
        for (idSeqDb,lengths) in self.seqLengthsAll():
            #Map IdSeqDbs from self to score ID either 
            #through remapping argument or through mapping
            #extracted from score meta-data argument
            if idSeqDb in idSeqDbToIdScoreRemapping:
                idScore = idSeqDbToIdScoreRemapping[idSeqDb]
            else:
                idScore = idSeqDbToIdScore[idSeqDb]
            #Map fragment ID through IdSeqDb mapping to score ID
            idFragToIdScore.update( 
                    dict(
                        ( (rec["id"],idScore) for rec in lengths )
                        )
                    )
        return idFragToIdScore


class ImmClassifierBenchMetricsSql(object):
    """Class that computes metrics for the taxonomic claasifier from pairwise test-predict SQL records"""

    ##input sample table
    sampTbl = "bench_samp"
    ##Data will be processed from this view name
    preFilterView = "bench_samp_flt"
    ##confusion table name
    confTbl = "conf"
    ##bottom-level taxa confusion table name
    confBotTbl = "conf_bot"
    ##model confusion table name
    confModTbl = "conf_mod"
    ##test count table name
    testCntTbl = "test_cnt"
    ##max test count table name
    testCntMaxTbl = "test_cnt_max"
    ##weight for the test count table name
    testCntWeightTbl = "test_cnt_wgt"
    ##weighted confusion table name
    confWeightedTbl = "conf_wgt"
    ##sensitivity table name with a group (i_lev_exc,i_lev_per,taxid_clade)
    sensTbl = "sens"
    ##specificity table name with a group (i_lev_exc,i_lev_per,taxid_clade)
    specTbl = "spec"
    ##sensitivity aggregate table name with a group (i_lev_exc,i_lev_per)
    sensAggrTbl = sensTbl+"_aggr"
    ##specificity aggregate table name with a group (i_lev_exc,i_lev_per)
    specAggrTbl = specTbl+"_aggr"
    ##accuracy aggregate table name with a group (i_lev_exc,i_lev_per)
    accuAggrTbl = "accu_aggr"
    ##suffix for euk tables
    eukSfx = "_euk"
    ##suffix for mic tables
    micSfx = "_mic"

    iLevNoExc = 0

    def __init__(self,db,taxaLevelsTbl,taxaNamesTbl):
        self.db = db
        self.taxaLevelsTbl = taxaLevelsTbl
        self.taxaNamesTbl = taxaNamesTbl

    def makeMetrics(self,nLevTestModelsMin=2,csvAggrOut="bench.csv",comment=None):
        self.preFilter(nLevTestModelsMin=nLevTestModelsMin)
        self.makeConfusion()
        self.makeConfusionBottomTaxa()
        self.makeConfusionModel()
        self.makeSampleCounts()
        self.makeSampleCountsWeight()
        self.makeConfusionWeighted()
        self.makeSensitivity()
        print "DEBUG: making specificity"
        self.makeSpecificity()
        print "DEBUG: making accuracy"
        self.makeAccuracy()
        self.makeAggrMetrics()
        self.repAggrMetrics(csvOut=csvAggrOut,comment=comment)

    def preFilter(self,nLevTestModelsMin):
        """Pre-filter data by creating a temporary view"""
        db = self.db
        db.ddl("""create temporary view %(preFilterView)s as 
                select * 
                from %(sampTbl)s
                where 
                    n_lev_test_models >= %(nLevTestModelsMin)s
                    or
                    i_lev_exc = %(iLevNoExc)s
                """ % dict(preFilterView=self.preFilterView,
                    nLevTestModelsMin=nLevTestModelsMin,
                    sampTbl=self.sampTbl,
                    iLevNoExc=self.iLevNoExc))

    def makeConfusionGeneric(self,sampTbl,confTbl,levTestFld,levPredFld,where=None):
        """Create a table with a confusion matrix in a sparse row_ind,col_ind,count format"""
        db = self.db
        if where is None:
            where = ""
        else:
            where = "where "+where
        db.createTableAs(confTbl,
                """select i_lev_exc,
                i_lev_per,
                %(levTestFld)s as id_lev_test,
                %(levPredFld)s as id_lev_pred,
                taxid_superking,
                count(*) as cnt
                from %(sampTbl)s 
                %(where)s
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    %(levTestFld)s,
                    %(levPredFld)s,
                    taxid_superking
                """ % dict(sampTbl=sampTbl,
                    levTestFld=levTestFld,
                    levPredFld=levPredFld,
                    where=where),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "id_lev_test",
                    "id_lev_pred",
                    "taxid_superking"]})
    
    def makeConfusion(self):
        """Create a table with a confusion matrix in a sparse row_ind,col_ind,count format"""
        self.makeConfusionGeneric(sampTbl=self.preFilterView,
                confTbl=self.confTbl,
                levTestFld="taxid_lev_test",
                levPredFld="taxid_lev_pred")

    def makeConfusionBottomTaxa(self):
        """Create a table with a confusion matrix for a bottom-most test and prediction taxa
        in a sparse row_ind,col_ind,count format"""
        self.makeConfusionGeneric(sampTbl=self.sampTbl,
                confTbl=self.confBotTbl,
                levTestFld="taxid_lev_test_bot",
                levPredFld="taxid_lev_pred_bot",
                where="i_lev_exc={i_lev} and i_lev_per={i_lev}".\
                        format(i_lev=self.iLevNoExc))

    def makeConfusionModel(self):
        """Create a table with a confusion matrix for a bottom-most test and prediction taxa
        in a sparse row_ind,col_ind,count format"""
        self.makeConfusionGeneric(sampTbl=self.sampTbl,
                confTbl=self.confModTbl,
                levTestFld="id_mod_test",
                levPredFld="id_mod_pred",
                where="i_lev_exc={i_lev} and i_lev_per={i_lev}".\
                        format(i_lev=self.iLevNoExc))

    def makeSampleCounts(self):
        """Create a table with sample counts per class"""
        db = self.db
        db.createTableAs(self.testCntTbl,
                """select i_lev_exc,
                i_lev_per,
                id_lev_test,
                sum(cnt) as sum_cnt
                from %(confTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    id_lev_test
                """ % dict(confTbl=self.confTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "id_lev_test"]})
    
    def makeSampleCountsWeight(self):
        """Create a table with weight that normalizes sample count by the max count in each group"""
        db = self.db
        
        ##find the max count within each group
        db.createTableAs(self.testCntMaxTbl,
                """select i_lev_exc,
                i_lev_per,
                max(sum_cnt) as max_sum_cnt
                from %(testCntTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per 
                """ % dict(testCntTbl=self.testCntTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per"]})
        
        ##compute the weight for each group as max_count/self-count
        db.createTableAs(self.testCntWeightTbl,
                """select 
                a.i_lev_exc as i_lev_exc,
                a.i_lev_per as i_lev_per,
                a.id_lev_test as id_lev_test,
                cast(b.max_sum_cnt as REAL)/a.sum_cnt as weight_cnt
                from 
                    %(testCntTbl)s a,
                    %(testCntMaxTbl)s b
                where
                    a.i_lev_exc = b.i_lev_exc
                    and
                    a.i_lev_per = b.i_lev_per
                """ % dict(testCntTbl=self.testCntTbl,
                    testCntMaxTbl=self.testCntMaxTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "id_lev_test"]})
    
    def makeConfusionWeighted(self):
        """Create a weighted confusion table.
        Specificity is sensitive to unbalanced
        testing set (small classes will be swamped
        by false positives from large classes).
        We balance the matrix by multiplying
        rows with weights that equalize the number
        of test cases across classes (as though we
        repeated the tests for small classes in
        order to match the larger ones)
        """
        db = self.db
        db.createTableAs(self.confWeightedTbl,
                """select a.i_lev_exc as i_lev_exc,
                a.i_lev_per as i_lev_per,
                a.id_lev_test as id_lev_test,
                a.id_lev_pred as id_lev_pred,
                a.taxid_superking as taxid_superking,
                a.cnt*b.weight_cnt as cnt
                from 
                    %(confTbl)s a,
                    %(testCntWeightTbl)s b
                where
                    a.i_lev_exc = b.i_lev_exc
                    and
                    a.i_lev_per = b.i_lev_per
                    and
                    a.id_lev_test = b.id_lev_test
                """ % dict(confTbl=self.confTbl,
                    testCntWeightTbl=self.testCntWeightTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "id_lev_test",
                    "id_lev_pred"]})
    
    def makeSensitivity(self):
        """Create a table with sensitivity metrics"""
        db = self.db
        db.createTableAs(self.sensTbl,
                """select i_lev_exc,
                i_lev_per,
                id_lev_test as id_lev,
                taxid_superking,
                sum((id_lev_test==id_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_clade
                from %(confTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    id_lev_test,
                    taxid_superking
                """ % dict(confTbl=self.confWeightedTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "id_lev",
                    "taxid_superking"]})

    def makeSpecificity(self):
        """Create a table with specificity metrics"""
        db = self.db
        #"where id_lev_pred > 0" excludes the reject group
        #"where id_lev_pred in (select ...)" only computes specificity
        #for clades that have positive testing samples, because otherwise we have
        #no chance at all to observe a specificity above zero.
        db.createTableAs(self.specTbl,
                """select i_lev_exc,
                i_lev_per,
                id_lev_pred as id_lev, 
                sum((id_lev_test=id_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_clade
                from %(confTbl)s a 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    id_lev_pred
                having 
                    id_lev_pred in 
                        (select 
                            id_lev_test 
                            from %(testCntTbl)s b
                            where 
                                a.i_lev_exc = b.i_lev_exc
                                and
                                a.i_lev_per = b.i_lev_per
                        )
                    and 
                    id_lev_pred > 0
                """ % dict(confTbl=self.confWeightedTbl,testCntTbl=self.testCntTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "id_lev"]})

    def makeAccuracy(self):
        """Create a table with accuracy metrics (over all samples).
        This is computed from non-weighted confusion matrix
        because on the weighted matrix it exactly equals sensitivity."""
        db = self.db
        db.createTableAs(self.accuAggrTbl,
                """select i_lev_exc,
                i_lev_per,
                sum((id_lev_test==id_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_lev
                from %(confTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per 
                """ % dict(confTbl=self.confTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per"]})
    
    def makeAggrMetrics(self):
        """Make tables with the aggregate metrics"""
        db = self.db
        for (tblPerCladeName,tblAggrName) in ((self.sensTbl,self.sensAggrTbl),
                (self.specTbl,self.specAggrTbl)):
            db.createTableAs(tblAggrName,
                    """select 
                    i_lev_exc,
                    i_lev_per,
                    avg(metr_clade) as metr_lev
                    from %(tblPerCladeName)s 
                    group by  
                        i_lev_exc,
                        i_lev_per 
                    """ % dict(tblPerCladeName=tblPerCladeName),
                    indices={"names":
                        ["i_lev_exc",
                        "i_lev_per"]})

    def repAggrMetrics(self,csvOut,comment=None):
        db = self.db
        taxaLevelsExtraTbl = "tmp_taxa_lev"
        db.createTableAs(taxaLevelsExtraTbl,
                """select id,level from %(taxaLevelsTbl)s
                """ % dict(taxaLevelsTbl=self.taxaLevelsTbl),
                indices={"names":["id"]})
        #db.ddl("""insert into %(taxaLevelsExtraTbl)s 
        #    (id,level) 
        #    values 
        #    (%(iLevNoExc)s,"none")
        #    """ % dict(taxaLevelsExtraTbl=taxaLevelsExtraTbl,
        #        iLevNoExc=self.iLevNoExc))
        for (iTbl,(tblAggrName,hdr)) in enumerate(
                (
                (self.sensAggrTbl,"Average per-clade sensitivity"),
                (self.specAggrTbl,"Average per-clade specificity"),
                (self.accuAggrTbl,"Per-sample accuracy")
                )
                ):
            sql = """
            select 
            b.level as lev_exc,
            c.level as lev_per,
            i_lev_exc,
            i_lev_per,
            a.metr_lev as metr_lev
            from 
                %(tblAggrName)s a,
                %(taxaLevelsTbl)s b,
                %(taxaLevelsTbl)s c
            where
                a.i_lev_exc = b.id
                and
                a.i_lev_per = c.id
            order by i_lev_exc,i_lev_per
            """ % dict(tblAggrName=tblAggrName,
                    taxaLevelsTbl=taxaLevelsExtraTbl)
            if iTbl == 0:
                mode = "w"
            else:
                mode = "a"
            if comment:
                fullComment = "%s (%s)" % (hdr,comment)
            else:
                fullComment = hdr
            db.exportAsPivotCsv(sql=sql,
                    out=csvOut,
                    mode=mode,
                    rowField="lev_exc",
                    colField="lev_per",
                    valField="metr_lev",
                    comment=fullComment,
                    colFieldOrderBy="i_lev_per",
                    restval="X",
                    valFormatStr="%.2f",
                    rowFieldOut="Exclude")

class ImmClassifierBenchToScore(object):
    
    sampTbl = ImmClassifierBenchMetricsSql.sampTbl

    def __init__(self):
        pass

    def catBenchSql(self,dbsInp,dbOut):
        if False:
            for (iInp,(lenSamp,dbInp)) in enumerate(sorted(dbsInp)):
                db = DbSqlLite(dbpath=dbInp,strategy="exclusive_unsafe")
                db.createIndices(self.sampTbl,names=["i_samp"])
                db.close()

        db = DbSqlLite(dbpath=dbOut,strategy="exclusive_unsafe")
        tableOutBase = "samp_acc"
        taxaLevels = TaxaLevels()
        for (iInp,(lenSamp,dbInp)) in enumerate(sorted(dbsInp)):
            db.ddl("attach database '%s' as dbin" % (dbInp,))
            #Using main. below is critical otherwise
            #db.creatTableAs wacks the dbin. table because
            #w/o the database scope SQLite resolves table name
            #into the last attached DB.
            tableOut = "main." + tableOutBase
            tableInp = "dbin." + self.sampTbl
            sqlArgs = dict(tableInp=tableInp,
                    lenSamp=lenSamp,
                    iLevExc=taxaLevels.getLevelId("species"),
                    iLevPerLowest=taxaLevels.getLevelId("genus"))
            sql = """select i_samp,
            %(lenSamp)s as len_samp,
            i_lev_per,
            score,
            taxid_lev_test=taxid_lev_pred as is_correct 
            from %(tableInp)s 
            where i_lev_exc=%(iLevExc)s 
            and i_samp in (
                select i_samp 
                from %(tableInp)s
                where i_lev_exc=%(iLevExc)s 
                    and i_lev_per=%(iLevPerLowest)s
                ) 
            order by i_samp,i_lev_per
            """ % sqlArgs
            if iInp == 0:
                db.createTableAs(tableOut,sql)
            else:
                db.ddl("""insert into %s
                %s""" % (tableOut,sql))

            db.ddl("detach database dbin")
        db.createIndices(tableOutBase,names=["i_samp","len_samp","i_lev_per","is_correct"])
        db.close()
