"""Benchmark support for ImmClassifierApp"""

from MGT.DirStore import *
from MGT.FastaIO import *
from Taxa import *

class ImmClassifierBenchmark(DirStore):

    fastaSfx = ".fna"
    
    def __init__(self,*l,**kw):
        DirStore.__init__(self,*l,**kw)
        self.seqDb = None
    
    def getFastaPath(self,idDb):
        return self.getFilePath("%s%s" % (idDb,self.fastaSfx))

    def listDbIds(self,iterPaths=None):
        """List DB IDs either from this store or from the externally provided iterable"""
        return list(self.fileNames(pattern="*"+self.fastaSfx,sfxStrip=self.fastaSfx,iterPaths=iterPaths))

    def getSeqDbIds(self):
        return set((str(x) for x in self.seqDb.getIdList()))
    
    def getImmIds(self,immDbs):
        immIds = set()
        for immDb in immDbs:
            immIds |= set((str(x) for x in immDb.listTaxids()))
        return immIds

    def selectIdsDb(self,immDbs=None):
        """Pick those IDs from SeqDb that will be used to build the benchmark.
        @param immDbs sequence of ImmDb instances - this is used to pick only those entries
        in SeqDb that also have models built against them. The self-models will
        still be excluded during testing - using ImmDb is just an easy way to
        filter out from SeqDb viruses as it was done for ImmDb. If None, no
        filtering by model IDs will be done."""
        ids = self.getSeqDbIds()
        if immDbs: 
            ids &= self.getImmIds(immDbs)
        return ids

    def makeSample(self,idDb,fragLen,fragCountMax):
        """Make a sample FASTA file for a given SeqDb ID"""
        outFasta = self.getFastaPath(idDb)
        self.shredFasta(idDb=idDb,outFasta=outFasta,fragLen=fragLen,
                fragCountMax=fragCountMax)

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
            inpFile = openCompressed(self.getFastaPath(idDb),"r")
            shutil.copyfileobj(inpFile,outFile,1024*1024)
            inpFile.close()
        if outClose:
            outFile.close()

    def shredFasta(self,idDb,outFasta,fragLen,
            fragCountMax,lineLen=80,outMode="w"):
        """Shred each record in multi-FASTA file into multiple records of fixed size"""
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
        outSeq = FastaWriter(out=outFasta,lineLen=lineLen,mode=outMode)
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
                idF = "%s_%s pos=%s nincat=%s len=%s" % (idDb,ind+iF,seqCatStart+stF,iF,len(sF))
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
                    i_lev_exc = 0
                """ % dict(preFilterView=self.preFilterView,
                    nLevTestModelsMin=nLevTestModelsMin,
                    sampTbl=self.sampTbl))

    def makeConfusionGeneric(self,sampTbl,confTbl,taxidLevTestFld,taxidLevPredFld):
        """Create a table with a confusion matrix in a sparse row_ind,col_ind,count format"""
        db = self.db
        db.createTableAs(confTbl,
                """select i_lev_exc,
                i_lev_per,
                %(taxidLevTestFld)s as taxid_lev_test,
                %(taxidLevPredFld)s as taxid_lev_pred,
                taxid_superking,
                count(*) as cnt
                from %(sampTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    %(taxidLevTestFld)s,
                    %(taxidLevPredFld)s,
                    taxid_superking
                """ % dict(sampTbl=sampTbl,
                    taxidLevTestFld=taxidLevTestFld,
                    taxidLevPredFld=taxidLevPredFld),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid_lev_test",
                    "taxid_lev_pred",
                    "taxid_superking"]})
    
    def makeConfusion(self):
        """Create a table with a confusion matrix in a sparse row_ind,col_ind,count format"""
        self.makeConfusionGeneric(sampTbl=self.preFilterView,
                confTbl=self.confTbl,
                taxidLevTestFld="taxid_lev_test",
                taxidLevPredFld="taxid_lev_pred")

    def makeConfusionBottomTaxa(self):
        """Create a table with a confusion matrix for a bottom-most test and prediction taxa
        in a sparse row_ind,col_ind,count format"""
        self.makeConfusionGeneric(sampTbl=self.sampTbl,
                confTbl=self.confBotTbl,
                taxidLevTestFld="taxid_lev_test_bot",
                taxidLevPredFld="taxid_lev_pred_bot")

    def makeSampleCounts(self):
        """Create a table with sample counts per class"""
        db = self.db
        db.createTableAs(self.testCntTbl,
                """select i_lev_exc,
                i_lev_per,
                taxid_lev_test,
                sum(cnt) as sum_cnt
                from %(confTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    taxid_lev_test
                """ % dict(confTbl=self.confTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid_lev_test"]})
    
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
                a.taxid_lev_test as taxid_lev_test,
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
                    "taxid_lev_test"]})
    
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
                a.taxid_lev_test as taxid_lev_test,
                a.taxid_lev_pred as taxid_lev_pred,
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
                    a.taxid_lev_test = b.taxid_lev_test
                """ % dict(confTbl=self.confTbl,
                    testCntWeightTbl=self.testCntWeightTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid_lev_test",
                    "taxid_lev_pred"]})
    
    def makeSensitivity(self):
        """Create a table with sensitivity metrics"""
        db = self.db
        db.createTableAs(self.sensTbl,
                """select i_lev_exc,
                i_lev_per,
                taxid_lev_test as taxid,
                taxid_superking,
                sum((taxid_lev_test==taxid_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_clade
                from %(confTbl)s 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    taxid_lev_test,
                    taxid_superking
                """ % dict(confTbl=self.confWeightedTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid",
                    "taxid_superking"]})

    def makeSpecificity(self):
        """Create a table with specificity metrics"""
        db = self.db
        #"where taxid_lev_pred > 0" excludes the reject group
        #"where taxid_lev_pred in (select ...)" only computes specificity
        #for clades that have positive testing samples, because otherwise we have
        #no chance at all to observe a specificity above zero.
        db.createTableAs(self.specTbl,
                """select i_lev_exc,
                i_lev_per,
                taxid_lev_pred as taxid, 
                sum((taxid_lev_test=taxid_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_clade
                from %(confTbl)s a 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    taxid_lev_pred
                having 
                    taxid_lev_pred in 
                        (select 
                            taxid_lev_test 
                            from %(testCntTbl)s b
                            where 
                                a.i_lev_exc = b.i_lev_exc
                                and
                                a.i_lev_per = b.i_lev_per
                        )
                    and 
                    taxid_lev_pred > 0
                """ % dict(confTbl=self.confWeightedTbl,testCntTbl=self.testCntTbl),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid"]})

    def makeAccuracy(self):
        """Create a table with accuracy metrics (over all samples).
        This is computed from non-weighted confusion matrix
        because on the weighted matrix it exactly equals sensitivity."""
        db = self.db
        db.createTableAs(self.accuAggrTbl,
                """select i_lev_exc,
                i_lev_per,
                sum((taxid_lev_test==taxid_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_lev
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
