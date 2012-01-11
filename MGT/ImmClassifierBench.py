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
    
    def selectIdsDb(self,immDbs):
        """Pick those IDs from SeqDb that will be used to build the benchmark.
        @param immDbs sequence of ImmDb instances - this is used to pick only those entries
        in SeqDb that also have models built against them. The self-models will
        still be excluded during testing - using ImmDb is just an easy way to
        filter out from SeqDb viruses as it was done for ImmDb"""
        immIds = set()
        for immDb in immDbs:
            immIds |= set((str(x) for x in immDb.listImmIds()))
        seqIds = set((str(x) for x in self.seqDb.getIdList()))
        return list(immIds & seqIds)

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
            #PreprocShredder treats sampNum()=0 as "all samples", so we need to set it
            #to at least 1, otherwise for catLen in a huge genome we will be getting
            #all samples. Maybe shredder's behaviour should be changed to sampNum<0 => all.
            #@todo Pick coords on a virtual concatenation of all sequences, otherwise
            #it will never be quite right.
            sampNum = lambda lab,seq,id: int(rndRound(len(seq)*seqLenRatio/fragLen))
        else:
            sampNum = 0
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

    ##confusion table name
    confTbl = "conf"
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

    def __init__(self,db,taxaLevelsTbl,taxaNamesTbl):
        self.db = db
        self.taxaLevelsTbl = taxaLevelsTbl
        self.taxaNamesTbl = taxaNamesTbl

    def makeMetrics(self,nLevTestModelsMin=2,csvAggrOut="bench.csv"):
        self.makeConfusion()
        self.makeSensitivity(nLevTestModelsMin=nLevTestModelsMin)
        self.makeSpecificity(nLevTestModelsMin=nLevTestModelsMin)
        self.makeAccuracy(nLevTestModelsMin=nLevTestModelsMin)
        self.makeAggrMetrics()
        self.repAggrMetrics(csvOut=csvAggrOut)

    def makeConfusion(self):
        """Create a table with a confusion matrix in a sparse row_ind,col_ind,count format"""
        db = self.db
        db.createTableAs(self.confTbl,
                """select i_lev_exc,
                i_lev_per,
                taxid_lev_test,
                taxid_lev_pred,
                taxid_superking,
                n_lev_test_models,
                count(*) as cnt
                from bench_samp 
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    taxid_lev_test,
                    taxid_lev_pred,
                    taxid_superking,
                    n_lev_test_models
                """,
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid_lev_test",
                    "taxid_lev_pred",
                    "taxid_superking",
                    "n_lev_test_models"]})
    
    def makeSensitivity(self,nLevTestModelsMin):
        """Create a table with sensitivity metrics"""
        db = self.db
        db.createTableAs(self.sensTbl,
                """select i_lev_exc,
                i_lev_per,
                taxid_lev_test as taxid,
                taxid_superking,
                sum((taxid_lev_test==taxid_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_clade
                from %(confTbl)s 
                where 
                    n_lev_test_models >= %(nLevTestModelsMin)s
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    taxid_lev_test,
                    taxid_superking
                """ % dict(confTbl=self.confTbl,nLevTestModelsMin=nLevTestModelsMin),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid",
                    "taxid_superking"]})

    def makeSpecificity(self,nLevTestModelsMin):
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
                from %(confTbl)s 
                where 
                    taxid_lev_pred in 
                        (select 
                            taxid_lev_test 
                            from %(confTbl)s
                            where 
                                n_lev_test_models >= %(nLevTestModelsMin)s
                        )
                    and 
                    taxid_lev_pred > 0
                group by  
                    i_lev_exc,
                    i_lev_per, 
                    taxid_lev_pred
                """ % dict(confTbl=self.confTbl,nLevTestModelsMin=nLevTestModelsMin),
                indices={"names":
                    ["i_lev_exc",
                    "i_lev_per",
                    "taxid"]})

    def makeAccuracy(self,nLevTestModelsMin):
        """Create a table with accuracy metrics (over all samples)"""
        db = self.db
        db.createTableAs(self.accuAggrTbl,
                """select i_lev_exc,
                i_lev_per,
                sum((taxid_lev_test==taxid_lev_pred)*cnt)/cast(sum(cnt) as REAL) as metr_lev
                from %(confTbl)s 
                where 
                    n_lev_test_models >= %(nLevTestModelsMin)s
                group by  
                    i_lev_exc,
                    i_lev_per 
                """ % dict(confTbl=self.confTbl,nLevTestModelsMin=nLevTestModelsMin),
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

    def repAggrMetrics(self,csvOut):
        db = self.db
        for (iTbl,(tblAggrName,comment)) in enumerate(
                (
                (self.sensAggrTbl,"Average per-clade sensitivity"),
                (self.specAggrTbl,"Average per-clade specificity"),
                (self.accuAggrTbl,"Accuracy")
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
                    taxaLevelsTbl=self.taxaLevelsTbl)
            if iTbl == 0:
                mode = "w"
            else:
                mode = "a"
            db.exportAsPivotCsv(sql=sql,
                    out=csvOut,
                    mode=mode,
                    rowField="lev_exc",
                    colField="lev_per",
                    valField="metr_lev",
                    comment=comment,
                    colFieldOrderBy="i_lev_per",
                    restval="X",
                    valFormatStr="%.2f")
            
