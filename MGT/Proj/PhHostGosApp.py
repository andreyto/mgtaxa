### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""GOS-specific prediction of bacterial hosts for the viral sequences.
"""

from MGT.App import *

from MGT.Taxa import *
from MGT.PhageHostApp import *
from MGT.SeqDbFasta import *

from MGT.Sql import *
from MGT.Svm import checkSaneAlphaHist

#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

class PhHostGosApp(App):
    """App-derived class GOS-specific prediction of bacterial hosts for the viral sequences"""

    batchDepModes = ("score-imms-gos","train-imms-gos")

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        topWorkDir = os.environ["GOSII_WORK"]
        opt.workDir = pjoin(topWorkDir,"ph-gos-bac")
        opt.cwd = opt.workDir
        opt.immDbRef = pjoin(topWorkDir,"icm-refseq")
        opt.immDbGos = pjoin(topWorkDir,"icm-gos-bac")
        opt.immIdToSeqIdsGos = pjoin(opt.immDbGos,"imm-seq-ids.pkl")
        opt.taxaDirGos = pjoin(topWorkDir,"ph-gos-bac-taxa")
        opt.taxaFileGos = pjoin(opt.taxaDirGos,"taxaTree.pkl")
        opt.seqDbGos = pjoin(opt.taxaDirGos,"seqdb")
        opt.scaffGos = "/usr/local/projects/GOSII/dougAnalysis/nonViral40kbPlusScaffolds.fa.gz"
        opt.scaffApis = "/usr/local/projects/GOSII/dougAnalysis/apisTaxonomicScaffoldAssignments.out.gz"
        opt.annotCont = '/usr/local/projects/GOSII/hlorenzi/IO/INTERACTIONS/contigs_functional_annotation.tab'
        #topPredDir = pjoin(topWorkDir,"ph-pred-random-inp-shuffle")
        #topPredDir = pjoin(topWorkDir,"ph-pred-random-inp-uni")
        topPredDir = pjoin(topWorkDir,"ph-pred")
        #opt.predSeq = "/usr/local/projects/GOSII/shannon/Indian_Ocean_Viral/asm_combined_454_large/454LargeContigs.fna"
        #opt.predSeq = pjoin(topWorkDir,"scaff-rnd","query.5K.fna")
        opt.predSeq = pjoin(topWorkDir,"scaff-gos-vir","asm_combined_454_large.5K.fna")
        opt.predOutDirRef = pjoin(topPredDir,"asm_combined_454_large")
        opt.predOutDirGos = pjoin(topPredDir,"asm_combined_454_large-gos-bac")
        opt.predOutDir = pjoin(topPredDir,"asm_combined_454_large-gos-bac-comb")
        opt.predOutMgtGosCsv = pjoin(topPredDir,"gos-bac_refseq","pred-taxa.csv")
        opt.predOutTaxaMetaCsv = pjoin(opt.predOutDir,"pred-taxa.meta.csv")
        opt.predOutDbSqlite = pjoin(opt.predOutDir,"pred.db.sqlite")
        opt.predOutStatsDir = pjoin(opt.predOutDir,"stats")
        opt.predOutStatsCsv = pjoin(opt.predOutStatsDir,"stats.csv")
        opt.nImmBatches = 200
        opt.predMinLenSamp = 5000
        opt.trainMinLenSamp = 100000
        opt.newTaxidTop = mgtTaxidFirst
        opt.newTaxNameTop = "jcvi_env_gos"
   
    def initWork(self,**kw):
        self.taxaTree = None
        pass
    
    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "make-custom-seq":
            return self.makeCustomTaxaTreeAndSeqDb(**kw)
        elif opt.mode == "train-imms-gos":
            return self.trainImmsGos(**kw)
        elif opt.mode == "score-imms-gos":
            return self.scoreImmsGos(**kw)
        elif opt.mode == "proc-scores":
            return self.processImmScores(**kw)
        elif opt.mode == "export-pred":
            return self.exportPred(**kw)
        elif opt.mode == "stats-pred":
            return self.statsPred(**kw)
        else:
            raise ValueError("Unknown mode: %s" % (opt.mode,))

    def makeCustomTaxaTreeAndSeqDb(self,**kw):
        """Add metagenomic sequences as custom nodes to the taxonomy tree and save the tree, and also create SeqDbFasta for them"""
        opt = self.opt
        rmrf(opt.taxaDirGos)
        makedir(opt.taxaDirGos)
        rmrf(opt.seqDbGos)
        makedir(opt.seqDbGos)
        seqDbGos = SeqDbFasta.open(path=opt.seqDbGos)
        compr = SymbolRunsCompressor(sym="N",minLen=1)
        nonDegenSymb = "ATCGatcg"
        gosTaxaTop = TaxaNode(id=opt.newTaxidTop,name=opt.newTaxNameTop,rank=unclassRank,divid=dividEnv,names=list())
        nextNewTaxid = gosTaxaTop.id + 1
        for rec in FastaReader(opt.scaffGos).records():
            hdr = rec.header()
            seqid = hdr.strip().split('>')[1] # scf768709870
            seq = compr(rec.sequence())
            if not checkSaneAlphaHist(seq,nonDegenSymb,minNonDegenRatio=0.99):
                print "WARNING: ratio of degenerate symbols is too high, "+\
                        "skipping the reference scaffold id %s" % (seqid,)
            if len(seq) >= opt.trainMinLenSamp:
                taxaNode = TaxaNode(id=nextNewTaxid,name=opt.newTaxNameTop+"_"+seqid,rank=unclassRank,divid=dividEnv,names=list())
                taxaNode.setParent(gosTaxaTop)
                fastaWriter = seqDbGos.fastaWriter(id=nextNewTaxid,lineLen=80)
                fastaWriter.record(header=hdr,sequence=seq)
                fastaWriter.close()
                nextNewTaxid += 1
        print "DEBUG: written %s sequence files in SeqDbFasta" % (nextNewTaxid - gosTaxaTop.id -1,)
        taxaTree = self.getTaxaTreeNcbi()
        envBacTop = taxaTree.getNode(envBacTaxid)
        gosTaxaTop.setParent(envBacTop)
        taxaTree.rebuild()
        taxaTreeStore = NodeStoragePickle(opt.taxaFileGos)
        taxaTreeStore.save(taxaTree)

    def trainImmsGos(self,**kw):
        """Train IMMs on bacterial scaffolds from the metagenome.
        For metagenomic scaffolds, we only train for leaf sequences. This
        is why we use ImmApp directly instead of higher level ImmClassifierApp"""
        opt = self.opt
        optI = copy(opt)
        optI.mode = "train"
        optI.seqDb = opt.seqDbGos
        optI.immDb = opt.immDbGos
        optI.immIdToSeqIds = opt.immIdToSeqIdsGos
        rmrf(opt.immDbGos)
        makedir(opt.immDbGos)
        dumpObj(makeDefaultImmSeqIds(optI.seqDb),optI.immIdToSeqIds)
        imm = ImmApp(opt=optI)
        return imm.run(**kw)

    def scoreImmsGos(self,**kw):
        """Score against IMMs trained on bacterial scaffolds from the metagenome.
        For metagenomic scaffolds, we only train for leaf sequences. This
        is why we use ImmApp directly instead of higher level ImmClassifierApp"""
        opt = self.opt
        optI = copy(opt)
        optI.mode = "score"
        optI.immDb = opt.immDbGos
        optI.immIdToSeqIds = opt.immIdToSeqIdsGos
        optI.outDir = opt.predOutDirGos
        optI.inpSeq = opt.predSeq
        rmrf(optI.outDir)
        makedir(optI.outDir)
        imm = ImmApp(opt=optI)
        return imm.run(**kw)
    
    @classmethod
    def _outScoreCombPath(klass,predOutDir):
        o = Struct(outDir=predOutDir)
        ImmClassifierApp.fillWithDefaultOptions(o)
        return o.outScoreComb
    
    def processImmScores(self,**kw):
        """Predict host taxonomy from IMM scores from RefSeq and Env contig models"""
        opt = self.opt
        makedir(opt.predOutDir)
        # merge scores from reference and metagenomic models
        scoreFileRef = self._outScoreCombPath(opt.predOutDirRef)
        scoreFileGos = self._outScoreCombPath(opt.predOutDirGos)
        scoreFileComb = self._outScoreCombPath(opt.predOutDir)
        scoreRef = loadObj(scoreFileRef)
        scoreGos = loadObj(scoreFileGos)
        # calling class method:
        scoreComb = scoreRef.catImms((scoreRef,scoreGos))
        dumpObj(scoreComb,scoreFileComb)
        jobs = self.predictByScore(scoreFileRef,scoreFileGos,scoreFileComb,**kw)

    def predictByScore(self,scoreFileRef,scoreFileGos,scoreFileComb,**kw):
        """Predict host by using a given score file and export predictions."""
        opt = self.opt
        # predict taxonomy using combined scores and customized tree
        scoreFiles = (scoreFileRef,scoreFileGos,scoreFileComb)
        predOutTaxaCsvFiles = [ scoreFile+'.csv' for scoreFile in scoreFiles ]
        jobs = []
        for (scoreFile,predOutTaxaCsvFile) in zip(scoreFiles,predOutTaxaCsvFiles):
            optI = copy(opt)
            optI.mode = "proc-scores"
            optI.outScoreComb = scoreFile
            optI.taxaTreePkl = opt.taxaFileGos
            optI.predOutTaxaCsv = predOutTaxaCsvFile
            PhageHostApp.fillWithDefaultOptions(optI)
            app = ImmClassifierApp(opt=optI)
            jobs+=app.run(**kw)
        jobs = []

        # load protein annotation and make annotation graphs
        ##optI.outScoreComb = scoreFileComb
        ##optI.mode = "cmp-prot-annot"
        ##app = PhageHostApp(opt=optI)

        ##jobs = app.run(depend=jobs)

        optJ = copy(opt)
        optJ.mode = "export-pred"
        optJ.predOutTaxaCsvRef, optJ.predOutTaxaCsvGos, optJ.predOutTaxaCsvComb = predOutTaxaCsvFiles
        app = self.factory(opt=optJ)

        jobs = app.run(depend=jobs)

        optJ.mode = "stats-pred"
        app = self.factory(opt=optJ)

        jobs = app.run(depend=jobs)
        return jobs

    def exportPred(self,**kw):
        """Export predictions for downstream analysis.
        This loads CSV files from PhageHostApp and joins it 
        with APIS assignments for host metagenomic sequences,
        and then saves as CSV again.
        Parameters are taken from self.opt.
        @param predOutTaxaCsvRef CSV prediction file produced by PhageHostApp based only on Ref DB
        @param predOutTaxaCsvGos CSV prediction file produced by PhageHostApp based only on Gos scaffolds
        @param predOutTaxaCsvComb CSV prediction file produced by PhageHostApp based on both Ref DB and Gos scaffolds
        @param predOutTaxametaCsv CSV prediction file for final output
        @param taxaFileGos custom taxonomy tree with metagenomic scaffolds/contigs
        """
        opt = self.opt
        self.loadMetaRefSql()
        self.loadMgtTaxaRefSql()
        db = DbSqlLite(dbpath=opt.predOutDbSqlite)
        csvOutFiles = (opt.predOutTaxaCsvRef, opt.predOutTaxaCsvGos, opt.predOutTaxaCsvComb)
        for (csvOutFile,tableSfx) in zip(csvOutFiles,("ref","met","comb")):
            db.createTableFromCsv(name="scaff_pred_"+tableSfx,
                    csvFile=csvOutFile,
                    hasHeader=True,
                    indices={"names":("name",)})
            db.createTableAs("pred_annot_"+tableSfx,
                    """select a.*,b.taxid as mgt_taxid,b.name as mgt_name,
                    b.name_genus as mgt_name_genus,
                    c.*
                    from scaff_pred_%s a
                    left outer join scaff_apis c
                    on a.name = c.id_seq
                    left outer join scaff_pred_mgt b
                    on a.name = b.id""" % (tableSfx,),
                    indices={"names":("id","name","taxid")})
            db.exportAsCsv("select * from pred_annot_"+tableSfx,opt.predOutTaxaMetaCsv+'_'+tableSfx)
        db.close()


    def loadMgtTaxaRefSql(self):
        """Load MGTAXA taxonomic assignments for metagenomic bacterial references into SQL tables.
        """
        opt = self.opt
        db = DbSqlLite(dbpath=opt.predOutDbSqlite)
        def _preProcFilter(row,fields,fieldInd):
            """Add the prefix to scaffold id"""
            id = "jcvi_env_gos_"+row[fieldInd["id"]]
            row[fieldInd["id"]] = id
            return (row,)
        db.createTableFromCsv(name="scaff_pred_mgt",
                csvFile=opt.predOutMgtGosCsv,
                hasHeader=True,
                preProc=_preProcFilter,
                fieldsMap={"len":SqlField(type="integer")},
                indices={"names":("id","taxid","name",)})
        db.close()

    def loadMetaRefSql(self):
        """Load annotations for metagenomic bacterial references into SQL tables.
        Currently loads file with APIS predictions projected onto scaffolds.
        """
        opt = self.opt
        db = DbSqlLite(dbpath=opt.predOutDbSqlite)
        apisInp = openCompressed(opt.scaffApis,"r")
        apisInp.next() # skip directory name at the top
        #scaffIds = set(db.selectAs1Col("select distinct name from scaff_pred"))
        #def _preProcFilter1(row,fields,fieldInd,scaffIds=scaffIds):
        #    """Add the prefix to scaffold id and return the row if it is in lookup set"""
        #    id_seq = "jcvi_env_gos_"+row[fieldInd["id_seq"]]
        #    if id_seq in scaffIds:
        #        row[fieldInd["id_seq"]] = id_seq
        #        return (row,)
        #    else:
        #        return tuple()
        def _preProcFilter(row,fields,fieldInd,minLen=opt.trainMinLenSamp):
            """Add the prefix to scaffold id and return the row if the sequence is long enough"""
            id_seq = "jcvi_env_gos_"+row[fieldInd["id_seq"]]
            if int(row[fieldInd["len"]]) >= minLen:
                row[fieldInd["id_seq"]] = id_seq
                return (row,)
            else:
                return tuple()
        db.createTableFromCsv(name="scaff_apis_1",
                csvFile=apisInp,
                fieldsMap={0:SqlField(name="id_seq"),
                    1:SqlField(name="len",type="integer"),
                    2:SqlField(name="genes",type="integer"),
                    16:SqlField(name="name_fam"),
                    17:SqlField(name="ratio_fam",type="real"),
                    18:SqlField(name="genes_fam",type="integer"),
                    19:SqlField(name="name_gen"),
                    20:SqlField(name="ratio_gen",type="real"),
                    21:SqlField(name="genes_gen",type="integer")},
                hasHeader=False,
                preProc=_preProcFilter)
        db.createTableAs("scaff_apis",
                """select *,
                (ratio_fam/100*genes_fam)/genes as ratio_fam_abs, 
                (ratio_gen/100*genes_gen)/genes as ratio_gen_abs
                from scaff_apis_1""",
                indices={"names":("id_seq","name_fam","name_gen")})
        apisInp.close()
        db.close()

    def statsPred(self,**kw):
        """Create aggregate tables and csv files to show various relationships between predictions"""
        from textwrap import wrap,dedent
        opt = self.opt
        db = DbSqlLite(dbpath=opt.predOutDbSqlite)
        makedir(opt.predOutStatsDir)
        outCsv = openCompressed(opt.predOutStatsCsv,'w')
        sqlAsComment = True
        epilog = "\n"
        sqlpar = dict(ratio_abs = 0.2)
        abstract = \
        """Abstract.
        Several tables that aggregate viral query sequences along various annotations of predicted
        bacterial hosts. Three tables are queried: pred_annot_ref where predictions were made against NCBI RefSeq bacteria
        (and archaea); pred_annot_met where predictions were made against sufficiently long metagenomic scaffolds classified
        as bacterial by external means (e.g. APIS gene votes); pred_annot_comb where predictions were made against a union
        of bacterial references used in the first two methods (in other words, NCBI + metagenomic scaffolds).
        The fields are named as follows: 
        mgt_* means metagenomic microbial references were classified by MGTAXA against bacterial NCBI sequences; mgt_name
        means whatever is the lowest taxonomic node for which the classification was made;
        name_fam, name_gen, etc are clades of the reference sequences predicted by APIS votes;
        ratio_fam_abs and ratio_gen_abs are ratios of the number of genes classified by APIS as name_fam and name_gen 
        respectively to the total number of genes on the metagenomic reference scaffold. We consider APIS assignment as
        "reliable" if either ratio is above %(ratio_abs)s.
        name,name_genus,name_species etc means taxonomic levels of NCBI reference sequences.\n"""
        outCsv.write(abstract)
        sql = """select mgt_name,count(*) as cnt from pred_annot_met group by mgt_name order by cnt desc""" % sqlpar
        comment = "Counts of phage classified against env reference grouped by MGTAXA name"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select mgt_name_genus,count(*) as cnt from pred_annot_met group by mgt_name_genus order by cnt desc""" % sqlpar
        comment = "Counts of phage classified against env reference grouped by MGTAXA genus name"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select a.*,b.* from pred_annot_ref a, pred_annot_comb b where a.id=b.id and not a.taxid=b.taxid and (b.ratio_gen_abs>=%(ratio_abs)s or b.ratio_fam_abs>=%(ratio_abs)s)""" % sqlpar
        db.exportAsCsv(sql,outCsv,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select name_gen,count(*) as cnt from scaff_apis where ratio_gen_abs>=%(ratio_abs)s group by name_gen order by cnt desc""" % sqlpar
        comment = "Counts of env reference grouped by APIS reliable genus (abundance of bacterial APIS clades, regardless of viral assignment)"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        # almost entirely to a single ref in gos
        sql = """select name_gen,count(*) as cnt from pred_annot_met where ratio_gen_abs>=%(ratio_abs)s group by name_gen order by cnt desc""" % sqlpar
        comment = "Counts of phage classified against env reference grouped by APIS reliable genus"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select name_family,count(*) as cnt from scaff_pred_ref group by name_family order by cnt desc"""
        comment = "Counts of phage classified against RefSeq bacteria grouped by family"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select name_genus,count(*) as cnt from scaff_pred_ref group by name_genus order by cnt desc"""
        comment = "Counts of phage classified against RefSeq bacteria grouped by genus"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select name_fam,count(*) as cnt from scaff_apis where ratio_fam_abs>=%(ratio_abs)s group by name_fam order by cnt desc""" % sqlpar
        comment = "Counts of env reference grouped by APIS reliable family"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select name_fam,count(*) as cnt from pred_annot_met where ratio_fam_abs>=%(ratio_abs)s group by name_fam order by cnt desc""" % sqlpar
        comment = "Counts of phage classified against env reference grouped by APIS reliable family"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select a.id,a.name_family,b.name_fam,b.id_seq from pred_annot_ref a, pred_annot_comb b where a.id=b.id and not a.taxid=b.taxid and b.ratio_fam_abs>=%(ratio_abs)s""" % sqlpar
        comment = "Query records for which env reference won over RefSeq reference and APIS family is reliable"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        sql = """select a.id,a.name_genus,b.name_gen,b.id_seq from pred_annot_ref a, pred_annot_comb b where a.id=b.id and not a.taxid=b.taxid and b.ratio_gen_abs>=%(ratio_abs)s""" % sqlpar
        comment = "Query records for which env reference won over RefSeq reference and APIS genus is reliable"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        # almost entirely to a single ref in gos
        sql = """select a.id,a.name_genus,b.name_gen,b.id_seq,b.ratio_gen_abs from pred_annot_ref a, pred_annot_met b where a.id=b.id and b.ratio_gen_abs>=%(ratio_abs)s order by b.name_gen,b.id_seq""" % sqlpar
        comment = "Query records for RefSeq and env reference assignments where APIS genus is reliable"
        db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog=epilog)
        outCsv.close()
        db.close()


    def getTaxaTreeNcbi(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(PhHostGosApp)
