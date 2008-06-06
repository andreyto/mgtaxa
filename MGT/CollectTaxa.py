"""Select from NCBI BLAST database the identifiers of sequence records with taxonomy assigned to them.
Create various summary tables with overview of the available data along the way.
"""

from MGT.Common import *
from MGT.Taxa import *
from MGT.Sql import *
from MGT.BlastDb import BlastDb
from SeqDb import *

class TaxaCollector(MGTOptions):
    """Collects data about taxonomically known sequence for training the classifier.
    Process the NCBI BLAST DB files by calling fastacmd.
    Aggregate, prioritize and partially remove redundancy along the taxonomy ids."""
    
    def __init__(self,dbSql):
        
        MGTOptions.__init__(self)
        self.dbSql = dbSql
        self.blastDb = BlastDb()
        self.taxaTree = None
        self.taxaLevels = None
        self.taxaTreeDbStore = NodeStorageDb(db=self.dbSql,tableSfx=self.taxaTreeTableSfxMain)
        

    def rebuild(self):
        #self.loadGiTaxNumpy()
        #self.loadRefseqAcc()
        #self.loadTaxTables()
        #self.loadSeq()
        #self.loadSeqToHdf()
        #self.indexHdfSeq()
        #self.selectActiveSeq()
        self.indexHdfActiveSeq()
        pass

    def loadTaxTables(self):
        """Create and fill all initial tables with taxa tree data.
        It only needs NCBU taxonomy dump files.
        It assigns dummy zero values to seq_len and seq_len_total attributes.
        Result: tree is saved in SQL DB"""
        self.loadTaxCategories()
        self.loadTaxNodes()
        self.loadTaxNames()
        self.loadTaxLevels()
        self.taxaTree.setAttribute("seq_len",0L)
        self.taxaTree.setAttribute("seq_len_tot",0L)
        self.taxaTreeDbStore.save(self.taxaTree)
        
    def loadSeq(self):
        db = self.dbSql
        self.loadGiTaxPickled()
        self.createTableSeq()
        #sql = """
                #load data infile '/usr/local/scratch/atovtchi.mgt_seq.bulk.tmp'
                #into table seq
                #fields
                #terminated by '|'
                #optionally enclosed by '"'
        #"""
        #sql = """
        #copy into seq from '/usr/local/scratch/atovtchi.mgt_seq.bulk.tmp'
        #"""
        #db.ddl(sql)

        self.idGenSeq = IntIdGenerator()
        inserterSeq        = db.makeBulkInserterFile(table='seq',bufLen=500000,workDir=self.tmpDir)
        inserterSeqMultiId = db.makeBulkInserterFile(table='seq_multi_id',bufLen=500000,workDir=self.tmpDir)
        inserterSeqHdr = db.makeBulkInserterFile(table='seq_hdr',bufLen=500000,workDir=self.tmpDir)
        blastDbs = self.blastDb.getDbs()
        self.loadSeqNCBI(blastDbs[0],inserterSeq,inserterSeqMultiId,inserterSeqHdr)
        self.loadSeqNCBI(blastDbs[1],inserterSeq,inserterSeqMultiId,inserterSeqHdr)
        self.loadSeqNCBI(blastDbs[2],inserterSeq,inserterSeqMultiId,inserterSeqHdr)
        self.loadSeqNCBI(blastDbs[3],inserterSeq,inserterSeqMultiId,inserterSeqHdr)
        self.loadSeqNCBI(blastDbs[4],inserterSeq,inserterSeqMultiId,inserterSeqHdr)
        inserterSeq.flush()
        inserterSeqMultiId.flush()
        inserterSeqHdr.flush()
        db.createIndices(table="seq",
        names=["gi","taxid","src_db","kind","project"],
        primary="id",
        compounds={"src":"taxid,src_db,kind,project"})
        db.ddl("analyze table seq",ifDialect="mysql")
        db.createIndices(table="seq_multi_id",
        names=["gi","acc_db"],
        primary="id_seq")
        db.ddl("analyze table seq_multi_id",ifDialect="mysql")
        db.createIndices(table="seq_hdr",
        names=["gi"],
        primary="id_seq")
        db.ddl("analyze table seq_hdr",ifDialect="mysql")
        self.clearGiTax()
        self.delDuplicateGiFromSeq()
        db.ddl("analyze table seq",ifDialect="mysql")
        db.ddl("analyze table seq_hdr",ifDialect="mysql")
        self.createFullTextIndexSeqHeader()

    def clearGiTax(self):
        self.gi2taxa = None
    
    def createTableSeq(self):
        # we removed "auto_increment primary key" from id
        # to speed up bulk load
        self.dbSql.ddl("""
        create table seq
        (
        id integer,
        gi bigint,
        taxid integer,
        src_db char(1),
        project char(4),
        seq_len bigint,
        acc_db char(3),
        acc char(18),
        kind char(2)
        )
        """,
        dropList=["table seq"])

        self.dbSql.ddl("""
        create table seq_hdr
        (
        id_seq integer,
        gi bigint,
        hdr varchar(%i)
        )
        """ % (self.fastaHdrSqlLen,),
        dropList=["table seq_hdr"])


        ## If we get a header record with two gi's like:
        ## >gi|23455713|ref|NC_004301.1| Enterobacteria phage FI, complete genome >gi|15183|emb|X07489.1| Bacteriophage SP genomic RNA
        ## then we insert the second gi along with db ("emb" here) into the seq_multi_id table, related to the record in 'seq' by
        ## "seq_multi_id.id_seq = seq.id".
        
        self.dbSql.ddl("""
        create table seq_multi_id
        (
        id_seq integer,
        gi bigint,
        acc_db char(3)
        )
        """,
        dropList=["table seq_multi_id"])
        
        
    def loadSeqNCBI(self,db,inserterSeq,inserterSeqMultiId,inserterSeqHdr):
        print "Processing BLAST DB " + db.db
        inp = self.blastDb.fastaReader(dbName=db.db,defLineTargetOnly=False)
        #curs = self.dbSql.cursor()
        iRec = 0
        #bufLen = 50000
        #sql = """
        #insert into seq
        #(id, gi, taxid, src_db, project, seq_len, acc_db, acc, kind, seq_hdr)
        #values
        #(%s,%s,%s,%s,%s,%s,%s,%s,%s)
        #"""
        #inserter = self.dbSql.makeBulkInserter(sql=sql,bufLen=bufLen)
        #inserter.n = 1
        #inserter.flush()
        #inp.close()
        #return
        for rec in inp.records():
            idSeq = self.idGenSeq()
            title = rec.header()[1:-1] #remove '>' and '\n'
            (gifld,gi,acc_db,acc,txt) = title.split('|',4)
            assert gifld == 'gi'
            gi = int(gi)
            defLine2 = txt.split(">gi|")
            if len(defLine2) > 1:
                txt = defLine2[0]
                (gi2,acc_db2,dummy) = defLine2[1].split('|',2)
                gi2 = int(gi2)
                inserterSeqMultiId((idSeq,gi2,acc_db2))
            try:
                taxid = int(self.gi2taxa[gi])
            except KeyError:
                print "Warning: Taxid not found for gi: "+gi
                taxid = 0
            # We are not interested in taxonomically unassigned sequence.
            # If we late change this, we also need to modify
            # exportIdsForSeqDb() where we currently do the inner join
            # with tax nodes table
            if taxid != 0:
                #Accession number formats are described here:
                #http://www.ncbi.nlm.nih.gov/Sequin/acc.html
                #WGS id might be 'AAAA00000000' or 'NZ_AAAA00000000', both cases seen in wgs db
                kind = ''
                acc_sfx = acc
                if acc_sfx[2] == '_':
                    kind = acc_sfx[:2]
                    acc_sfx = acc_sfx[3:]
                project = ''
                if acc_sfx[:4].isalpha() and acc_sfx[5].isdigit():
                    project = acc_sfx[:4]
                seqLen = rec.seqLen()
                values = (idSeq,gi,taxid,db.id,project,seqLen,acc_db,acc,kind)
                #values = [ str(x) for x in values ]
                inserterSeq(values)
                values = (idSeq,gi,title[:self.fastaHdrSqlLen])
                inserterSeqHdr(values)
            if iRec % 50000 == 0:
                print db.db, title, taxid, iRec, seqLen
                #if iRec >= 500000:
                #    break
            iRec += 1
        inp.close()

    def tmp_loadSeqHdr(self):
        idGenSeq = IntIdGenerator()
        db = self.dbSql
        db.ddl("""
        create table seq_hdr
        (
        id_seq integer,
        gi bigint,
        hdr varchar(%i)
        )
        """ % (self.fastaHdrSqlLen,),
        dropList=["table seq_hdr"])
        reCutLen = re.compile(r"(.*)\Wlen:[0-9]+")
        inserterSeqHdr = db.makeBulkInserterFile(table='seq_hdr',bufLen=500000,workDir=self.tmpDir)
        inp = open("all.hdr",'r')
        for rec in inp:
            idSeq = idGenSeq()
            title = reCutLen.search(rec[1:-1]).group(1) #remove '>' and '\n', then cut the " len:0000000" suffix
            (gifld,gi,txt) = title.split('|',2)
            assert gifld == 'gi'
            gi = int(gi)
            values = (idSeq,gi,title)
            inserterSeqHdr(values)
        inserterSeqHdr.flush()
        db.createIndices(table="seq_hdr",
        names=["gi"],
        primary="id_seq")
        db.ddl("analyze table seq_hdr",ifDialect="mysql")

    def createFullTextIndexSeqHeader(self):
        """Build a full text index for sequence FASTA hdeaders in 'seq_hdr' table (currently implemented only for MySQL back-end).
        Warning: it took 1 hr for 28M records."""
        db = self.dbSql
        db.ddl("alter table seq_hdr add fulltext index hdr(hdr)",ifDialect="mysql",dropList=["index hdr on seq_hdr"])
        

    def loadRefseqAcc(self):
        self.dbSql.ddl("""\
        create table refseq_acc
        (
        prefix char(2),
        acc char(40),
        molecule char(10),
        method char(15),
        descr varchar(400)
        )
        """,
        dropList=["table refseq_acc"])
        self.dbSql.createIndices(names=["prefix"],table="refseq_acc")
        
        self.dbSql.executemany("""\
        insert into refseq_acc
        (prefix,acc,molecule,method,descr)
        values
        (%s,%s,%s,%s,%s)
        """,refseqAccFormat)
        
    def loadTaxCategories(self):
        self.dbSql.ddl("""\
        create table taxa_cat
        (cat char(1), taxid_ancestor integer, taxid integer)
        """,
        dropList=["table taxa_cat"])
        inp = open(self.taxaCatFile,'r')
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 500000
        sql = """\
        insert into taxa_cat
        (cat,taxid_ancestor,taxid)
        values
        (%s,%s,%s)
        """
        inserter = self.dbSql.makeBulkInserterFile(table="taxa_cat",bufLen=500000,workDir=self.tmpDir)
        for rec in inp:
            vals = rec.split()
            vals[1] = int(vals[1])
            vals[2] = int(vals[2])
            inserter(vals)
        inp.close()
        inserter.flush()
        curs.close()
        self.dbSql.createIndices(names=["cat","taxid_ancestor"],table="taxa_cat",primary="taxid")

    def loadTaxNodes(self):
        self.dbSql.ddl("""\
        create table taxa_node
        (
        taxid integer,
        partaxid integer,
        rank char(20),
        embl_code  char(2),
        divid  integer,
        inh_div  bool,
        gcode_id integer,
        inh_gc  bool,
        mgcode_id integer,
        inhmgc  bool,
        gbhidden bool,
        hidsubtree bool,
        comments  char(40)
        )
        """,
        dropList=["table taxa_node"])
        inp = open(self.taxaNodesFile,'r')
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 500000
        sql = """\
        insert into taxa_node
        (
        taxid,
        partaxid,
        rank,
        embl_code,
        divid,
        inh_div,
        gcode_id,
        inh_gc,
        mgcode_id,
        inhmgc,
        gbhidden,
        hidsubtree,
        comments
        ) 
        values
        (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
        """
        inserter = self.dbSql.makeBulkInserterFile(table="taxa_node",bufLen=500000,workDir=self.tmpDir)
        for rec in inp:
            # The last '|' is behind the last field, so we cut it
            inserter([ x.strip() for x in rec.split('|')][:-1])
        inserter.flush()
        inp.close()
        curs.close()
        #Get rid of spaces in 'rank' field to simplify export of data later on
        self.dbSql.execute("update taxa_node set rank = replace(rank,' ','_')").close()
        self.dbSql.createIndices(names=["partaxid","divid"],table="taxa_node",primary="taxid")


    def loadTaxNames(self):
        """We load only 'scientific name' entries."""
        self.dbSql.ddl("""\
        create table taxa_names
        (
        taxid integer,
        name  varchar(180)
        )
        """,
        dropList=["table taxa_names"])
        inp = open(self.taxaNamesFile,'r')
        inserter = self.dbSql.makeBulkInserterFile(table="taxa_names",bufLen=500000,workDir=self.tmpDir)
        for line in inp:
            rec = [ x.strip() for x in line.split('|') ]
            if rec[3] == "scientific name":
                inserter((int(rec[0]),rec[1]))
        inserter.flush()
        inp.close()
        self.dbSql.createIndices(table="taxa_names",primary="taxid")
        self.dbSql.executeAndAssertEmpty("select taxid from taxa_node where taxid not in (select taxid from taxa_names)")

    def loadTaxLevels(self):
        self.loadTaxLevelsMem(withTaxTree=True)
        #self.taxaLevels.loadTaxLevelsRows()
        #self.taxaLevels.loadTaxLevelsColumns()
        #self.taxaLevels.makeStatsTables()
        
    def loadTaxLevelsMem(self,withTaxTree=False):
        if self.taxaLevels is None:
            if withTaxTree:
                self.loadTaxNodesMem()
            # that will also set level attributes inside tree nodes
            self.taxaLevels = TaxaLevelsDb(self.dbSql,self.taxaTree)
        
    def loadTaxNodesMem(self):
        if self.taxaTree is None:
            self.taxaTree = \
                TaxaTree(NodeStorageNcbiDump(ncbiDumpFile=self.taxaNodesFile,
                                             ncbiNamesDumpFile=self.taxaNamesFile))
        #self.taxaTree.write(sys.stdout)

       
    def loadGiTaxSql(self):
        self.dbSql.dropTable("gi_taxa")
        self.dbSql.execute(
        """
        create table gi_taxa
        (
        gi bigint,
        taxid integer
        )
        """)
        #self.dbSql.executemany("insert into gi_tax (gi,taxid) values (%s,%s)",[(1,2),(3,4)])
        #print "Test done"
        #self.dbSql.execute("insert into gi_tax (gi,taxid) values ('1','2')")
        inp = file(self.taxaGiFile,'r')
        iRow = 0
        bufLen = 100000
        sql = "insert into gi_taxa (gi,taxid) values (%s,%s)"
        curs = self.dbSql.cursor()
        inserter = self.dbSql.makeBulkInserterFile(table='gi_taxa',bufLen=500000,workDir=self.tmpDir)
        print "Start insert"
        for record in inp:
            inserter(record.split())
            iRow += 1
            if iRow % 100000 == 0:
                print record,iRow
        inserter.flush()
        curs.close()
        inp.close()
        print "Start index build"
        self.dbSql.ddl("""\
        alter table gi_taxa  
        ADD PRIMARY KEY gi(gi), 
        add index taxid(taxid)
        """)
        

    def loadGiTaxBdb(self,inpFile):
        import anydbm, whichdb
        print whichdb.whichdb('/export/atovtchi/taxa.db')
        gi2taxa = anydbm.open('/export/atovtchi/taxa.db', 'c')
        inp = file(inpFile,'r')
        iRow = 0
        buff = {}
        for line in inp:
            (gi,taxid) = line.split()
            buff[gi] = taxid
            if iRow % 100000 == 0:
                print gi, taxid, iRow
                gi2taxa.update(buff)
                buff = {}
            iRow += 1
        inp.close()
        taxaCnt = {}
        for (gi,taxid) in gi2taxa.iteritems():
            taxaCnt[taxid] = taxaCnt.get(taxid,0) + 1
        print sorted(((cnt,taxid) for (taxid,cnt) in taxaCnt.iteritems()))

    def loadGiTaxNumpy(self):
        inp = openCompressed(self.taxaGiFile,'r')
        #inp = open("test_gi_taxid_nucl.dmp",'r')
        #twocol = numpy.loadtxt(inp,dtype=numpy.int32)
        twocol = numpy.fromfile(inp,dtype="i4",sep='\n')
        twocol.shape = (twocol.shape[0]/2,2)
        gi2taxa = fromWhereItems({'ind':twocol[:,0], 'val':twocol[:,1]})
        print gi2taxa.shape, gi2taxa.dtype
        dumpObj(gi2taxa,self.taxaPickled)
        #pdb.set_trace() 


    def loadGiTaxPickled(self):
        self.gi2taxa = loadObj(self.taxaPickled)
        #pdb.set_trace()


    def delDuplicateGiFromSeq(self):
        """Delete duplicate records (by gi) from seq table, leave only those with smallest gi.
        Duplicate records appear because some sequences (with identical defline) are included
        into both 'refseq_genomic' and 'other_genomic' BLAST databases."""
        db = self.dbSql
        db.createTableAs("tmp_gi2","""
            (select gi,count(*) as cnt from seq group by gi having count(*) > 1)
        """)
        db.createIndices(names=["gi"],table="tmp_gi2")
        db.createTableAs("tmp_id_gi2","""
            (select a.id,a.gi from seq a,tmp_gi2 b where a.gi = b.gi)
        """)
        db.createIndices(names=["id","gi"],table="tmp_id_gi2")
        db.createTableAs("tmp_id_gi1","""
        (select a.id from tmp_id_gi2 a where a.id > SOME ( select b.id from tmp_id_gi2 b where a.gi = b.gi ))
        """)
        db.createIndices(names=["id"],table="tmp_id_gi1")
        db.ddl("""
            delete from seq where id in (select id from tmp_id_gi1)
        """)
        db.ddl("""
            delete from seq_hdr where id_seq in (select id from tmp_id_gi1)
        """)
        db.dropTables(("tmp_gi2","tmp_id_gi2","tmp_id_gi1"))

        db.createIndices(table="seq",
                         names=["gi"],
                         attrib={"gi":{"unique":True}})
        db.createIndices(table="seq_hdr",
                         names=["gi"],
                         attrib={"gi":{"unique":True}})

    def selectActiveSeq(self):
        """Select 'active' source sequence set - the one that can be used for training.
        @post tables act_seq and act_src contain usable training sequence."""
        self.excludePreSource()
        self.selectSeqBySource()
        self.excludePostSource()

    def excludePreSource(self):
        """Exclude some sequence before selectTaxSource()
        @post table seq_excl has excluded idseq's."""
        ## Currently we use MySQL specific full text search
        db = self.dbSql
        db.createTableAs("seq_excl","""
            (select id_seq as id from seq_hdr where MATCH(hdr) AGAINST('"ribosomal RNA" "rRNA"' IN BOOLEAN MODE) )
            """)
        db.createIndices(primary="id",table="seq_excl")

        
    def selectSeqBySource(self):
        """For each taxid, group available sequence by its source (refseq, wgs, htgs, nt) and select groups by priority.
        @post table all_src that has all selected groups and synthetic primary key 'id' for each group."""
        db = self.dbSql
        db.createTableAs("taxa_src_1","""
            (SELECT taxid                ,
                            src_db       ,
                            kind         ,
                            project      ,
                            0            AS priority,
                            COUNT(*)     AS cnt     ,
                            SUM(seq_len) AS seq_len
                    FROM    seq
                    WHERE   id not in (SELECT id FROM seq_excl)
                    GROUP BY taxid ,
                            src_db ,
                            kind   ,
                            project
            )
        """)

        db.createTableAs("taxa_src","""
                (SELECT *
                        FROM    taxa_src_1
                )
        ORDER BY taxid ,
                src_db ,
                kind   ,
                project
        """)
        
        db.dropTable("taxa_src_1")
        
        #db.ddl("""ALTER TABLE taxa_src ADD id INTEGER auto_increment PRIMARY KEY""",ifDialect="mysql")
        db.ddl("""ALTER TABLE taxa_src ADD id INTEGER auto_increment  PRIMARY KEY""")
        db.createIndices(names=["taxid","src_db","kind","project"],table="taxa_src")
        
        ## Drop "Alternative" full genomes (e.g. haplotypes)
        db.ddl("""
        DELETE
        FROM    taxa_src
        WHERE   kind = 'AC'
        """)

        db.createTableAs("wgs_src_1","""
        (SELECT *
                FROM    taxa_src
                WHERE   project <> ''
                    AND kind    <> ''
        )
        """)
        
        db.ddl("""
        INSERT
        INTO    wgs_src_1
        SELECT  a.*
        FROM    taxa_src a
        WHERE   a.project   <> ''
            AND a.taxid NOT IN
                (SELECT taxid
                FROM    wgs_src_1
                )
        """)

        db.createTableAs("wgs_src_2","""
        (SELECT a.*
                FROM    wgs_src_1 a
                WHERE   a.seq_len >= ALL
                        (SELECT b.seq_len
                        FROM    wgs_src_1 b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)

        db.createTableAs("wgs_src_3","""
        (SELECT a.*
                FROM    wgs_src_2 a
                WHERE   a.id >=
                        (SELECT MAX(b.id)
                        FROM    wgs_src_2 b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)
        db.dropTable("wgs_src_2")

        db.createTableAs("wgs_src_nr","""
        (SELECT *
                FROM    wgs_src_3
        )
        """)

        db.createTableAs("nc_src_1","""
        (SELECT a.*
                FROM    taxa_src a
                WHERE   a.kind        = 'NC'
                    AND a.seq_len * 2 > ALL
                        (SELECT b.seq_len
                        FROM    wgs_src_nr b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)
        
        db.createTableAs("htg_over_nc_src_1","""
        (SELECT a.*
                FROM    taxa_src a,
                        nc_src_1 b
                WHERE   a.src_db     = 'h'
                    AND a.kind      <> 'NC'
                    AND a.taxid      = b.taxid
                    AND a.seq_len    > b.seq_len * 5
                    AND b.seq_len    < 500000
                    AND a.taxid NOT IN
                        (SELECT taxid
                        FROM    wgs_src_nr
                        )
        )
        """)

        db.createTableAs("nc_src_2","""
        (SELECT *
                FROM    nc_src_1
                WHERE   taxid NOT IN
                        (SELECT taxid
                        FROM    htg_over_nc_src_1
                        )
        )
        """)

        db.createIndices(names=["taxid"],table="nc_src_2")
        
        db.createTableAs("nc_src_3","""
        (SELECT a.*
                FROM    nc_src_2 a
                WHERE   a.id >= 
                        (SELECT MAX(b.id)
                        FROM    nc_src_2 b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)
        
        db.createTableAs("nc_src","""
        (SELECT a.*,
                        'nc' AS stage
                FROM    nc_src_3 a
        )
        """)

        db.createIndices(names=["taxid"],table="nc_src_3")
        
        db.createTableAs("wgs_src","""
        (SELECT a.*,
                        'wg' AS stage
                FROM    wgs_src_nr a
                WHERE   taxid NOT IN
                        (SELECT taxid
                        FROM    nc_src
                        )
        )
        """)
        
        db.createIndices(names=["taxid"],table="wgs_src")
        
        db.createTableAs("htg_src_1","""
        (SELECT *
                FROM    taxa_src
                WHERE ( src_db      = 'h'
                     OR kind       <> '' )
                    AND (taxid NOT IN
                        (SELECT taxid
                        FROM    nc_src
                        )
                    AND taxid NOT IN
                        (SELECT taxid
                        FROM    wgs_src
                        ) )
        )
        """)
        
        db.createTableAs("htg_src","""
        (SELECT a.*,
                        'ht' AS stage
                FROM    htg_src_1 a
        )
        """)
        
        db.createTableAs("gen_src","""
        (SELECT *
                FROM    nc_src
        )
        """)
        
        db.ddl("""
        INSERT
        INTO    gen_src
        SELECT  *
        FROM    wgs_src
        
        UNION
        
        SELECT  *
        FROM    htg_src
        """)

        db.createIndices(names=["taxid"],table="gen_src")
        
        db.createTableAs("nt_src","""
        (SELECT a.*,
                        'nt' AS stage
                FROM    taxa_src a
                WHERE   taxid NOT IN
                        (SELECT taxid
                        FROM    gen_src
                        )
        )
        """)

        db.createTableAs("all_src_1","""
        (SELECT a.*,
                        'g' AS src_type
                FROM    gen_src a
        )
        """)
        
        db.ddl("""
        INSERT
        INTO    all_src_1
        SELECT  a.*,
                'o' AS src_type
        FROM    nt_src a
        """)
    
    #def selectTaxSet(self):

        #db = self.dbSql
        ## We exclude records that we do not need here.
        db.createTableAs("all_src","""
        (SELECT a.*            ,
                        b.cat  ,
                        c.divid,
                        c.rank
                FROM    all_src_1 a
                        LEFT JOIN taxa_cat b
                        ON      a.taxid = b.taxid
                        LEFT JOIN taxa_node c
                        ON      a.taxid = c.taxid
                WHERE   not (
                    a.taxid IS NULL
                    OR a.taxid = 0
                    OR c.divid IS NULL
                    OR b.cat   IS NULL
                    OR c.divid IN (7,8,11)
                    )
        )
        """)

        db.ddl("ANALYZE TABLE all_src",ifDialect="mysql")

        db.createIndices(names=["taxid","src_db","kind","project","cat","divid","rank"],
            table="all_src",
            primary="id",
            compounds={"src":"taxid,src_db,kind,project"},
            attrib={"src":{"unique":True}})
        
        db.ddl("ANALYZE TABLE all_src",ifDialect="mysql")


    def excludePostSource(self):
        """Exclude some sequence after selectTaxSource().
        @post table act_seq is mapping into seq and act_src
        @post table act_src is all_src with some records
        dropped and seq_len recomputed from seq_src.
        The compination of these two new tables allows to
        filter initial sequence set described by all_src 
        both at individual sequence level (e.g. drop short sequences)
        and at source id level (meaning origin,sequence type)
        (e.g. retain only longest RefSeq strain for each species).
        Currently this method does not drop any sequence.
        After this method, act_seq and act_src can be used
        to randomly sample for training and testing sets."""
        db = self.dbSql

        db.createTableAs("act_seq","""
        (SELECT         a.id,
                        b.id as id_src,
                        a.seq_len
                FROM    seq a,
                        all_src b
                WHERE   a.taxid   = b.taxid
                    AND a.src_db  = b.src_db
                    AND a.kind    = b.kind
                    AND a.project = b.project
                    AND a.taxid  <> 0
                    AND a.id not in (SELECT id FROM seq_excl)
        )
        """)
        db.createIndices(names=["id_src"],primary="id",table="act_seq")
        db.ddl("ANALYZE TABLE act_seq",ifDialect="mysql")

        db.createTableAs("act_src","""
        (SELECT a.id,a.taxid,a.src_db,a.kind,a.project,a.cat,a.divid,
                b.seq_len
                FROM    all_src a,
                        (SELECT id_src, sum(seq_len) as seq_len FROM act_seq GROUP BY id_src) b
                WHERE  a.id = b.id_src
        )
        """)

        db.ddl("ANALYZE TABLE act_src",ifDialect="mysql")

        db.createIndices(names=["taxid","src_db","kind","project","cat","divid"],
            table="act_src",
            primary="id",
            compounds={"src":"taxid,src_db,kind,project"},
            attrib={"src":{"unique":True}})
        
        db.ddl("ANALYZE TABLE act_src",ifDialect="mysql")

        # This is more for reporting - operational counts are based
        # on samples and created in LabelTaxa module
        db.createTableAs("taxa_seq_len","""
        select taxid,sum(seq_len) as seq_len from act_src group by taxid
        """)
        db.createIndices(primary="taxid",table="taxa_seq_len")
        db.ddl("ANALYZE TABLE taxa_seq_len",ifDialect="mysql")


    def exportIdsForSeqDb(self):
        """Write text files with sequence ids that will be used to create sequence HDF5 files.
        Currently we load only NCBI sequence, but we want to decouple SeqDB from NCBI GIs
        in case we will also use non-NCBI sequence in the future. Therefore, we write
        a file with (GI,ID) pairs where ID is our internal sequence id. The order is defined by
        our tree nested set index. That most closely corresponds to the order in which
        we will be traversing the tree most of the time."""
        treeTable = self.taxaTreeDbStore.tblNodes
        db = self.dbSql
        # GI,ID list for fastacmd input and HDF index id
        db.exportToFile(\
        """SELECT  a.gi,a.id""",
        """
        FROM            seq a,
                        %(treeTable)s b
        WHERE           a.taxid = b.id
        ORDER BY        b.lnest,a.id
        """ % {"treeTable":treeTable},
        fileName=self.seqGiIdFile,fieldsTerm=' ',linesTerm=r'\n')


    def selectSeqIds(self):

        db = self.dbSql

        self.loadTaxLevelsMem()

        db.createTableAs("seq_sel","""
        (SELECT         a.gi                                           ,
                        a.taxid                                        ,
                        b.id as src_id
                FROM    seq a,
                        all_src b
                WHERE   a.taxid   = b.taxid
                    AND a.src_db  = b.src_db
                    AND a.kind    = b.kind
                    AND a.project = b.project
                    AND a.taxid  <> 0
        )
        """)

        ## Proved to be essential in MySQL for efficient planning of queries
        db.ddl("ANALYZE TABLE seq_sel",ifDialect="mysql")

        db.createIndices(names=["taxid","src_id"],
            table="seq_sel",
            primary="gi")

        db.ddl("ANALYZE TABLE seq_sel",ifDialect="mysql")
        
        ### in MySQL computation of acc_no_ver (ACC w/o '.X' suffix)
        ### would be just SUBSTRING_INDEX( a.acc , '.', 1 ) AS acc_no_ver,
        ### but we use SQL standard conforming expression
        
        #db.createTableAs("seq_sel","""
        #(SELECT a.*                                                    ,
                        #b.cat                                          ,
                        #b.stage                                        ,
                        #b.src_type                                     ,
                        #SUBSTRING(a.acc FROM 1 FOR
                        #CASE
                           #WHEN POSITION('.' IN a.acc) <> 0
                           #THEN POSITION('.' IN a.acc) - 1
                           #ELSE LENGTH(a.acc)
                        #END)
                        #AS acc_no_ver                                  ,
                        #b.divid                                        ,
                        #b.rank
                #FROM    seq a,
                        #all_src b
                #WHERE   a.taxid   = b.taxid
                    #AND a.src_db  = b.src_db
                    #AND a.kind    = b.kind
                    #AND a.project = b.project
                    #AND a.taxid  <> 0
        #)
        #""")
        
        #db.createIndices(names=["taxid","src_db","kind","project","cat",
            #"stage","src_type","gi","acc","acc_no_ver","divid","rank"],
            #table="seq_sel",
            #primary="id")

        levelsComma = self.taxaLevels.getLevelColumnsComma(alias='b',order='descend')

        ## save gis only into a file to be used as gi list for NCBI alias db
        db.exportToFile(\
        """SELECT  a.gi""",
        """FROM    seq_sel a, taxa_level_col b
        WHERE a.taxid = b.taxid
        ORDER BY %s,a.taxid,a.gi
        """ % (levelsComma,),
        fileName=self.collectTaxaGiFile,fieldsTerm=' ',linesTerm=r'\n')

        ### save all important fields for selected seq records into a
        ### flat file to be merged with
        ### fasta sequence data extracted from blast databases
        ### through just created gi list
        ### the order must match the preceding 'select' for gis
        ### Full csv format will have something like OPTIONALLY ENCLOSED BY '"'

        #db.ddl("""
        #SELECT  gi      ,
                #taxid   ,
                #src_db  ,
                #kind    ,
                #project ,
                #cat     ,
                #stage   ,
                #src_type,
                #id     ,
                #seq_len ,
                #divid   ,
                #rank
        #INTO    OUTFILE '/home/atovtchi/scratch/mgtdata/phyla_sel.csv'
        #FIELDS TERMINATED BY ' '
        #LINES TERMINATED BY '\n'
        #FROM    seq_sel
        #ORDER BY taxid,
                #id
        #""")


    def loadSeqToHdf(self):
        """Load sequence data for all records in 'seq' table into HDF dataset."""
        self.exportIdsForSeqDb()
        ## CAUTION! mode="w" trancates the HDF file if it already exists.
        ## "a" opens file for writing w/o destroying existing data,
        ## "r+" is the same as "a" but file must exist
        seq_len_tot = long(self.dbSql.selectScalar("select sum(seq_len) from seq"))
        seq_cnt = long(self.dbSql.selectScalar("select count(*) from seq"))
        hdfFile = pt.openFile(self.hdfSeqFile,mode="w")
        hdfLoader = HdfSeqLoader(hdfFile=hdfFile,
                hdfGroup=self.hdfSeqGroup,
                seqSizeEstimate=seq_len_tot,
                seqNumEstimate=seq_cnt)
        hdfLoader.loadBlastDb(giIdFile=self.seqGiIdFile)
        #hdfLoader.close()
        hdfFile.close()

    def indexHdfSeq(self):
        """Create SQL tables that index sequence in HDF dataset."""
        db = self.dbSql
        hdfFile = pt.openFile(self.hdfSeqFile,mode="r")
        hdfSeq = HdfSeqReaderSql(db=db,hdfFile=hdfFile,hdfGroup=self.hdfSeqGroup)
        hdfSeq.loadIndToSql()
        #hdfSeq.close()
        hdfFile.close()

    def indexHdfActiveSeq(self):
        """Create HDF index dataset for active sequence (from act_seq table)."""
        hdfFile = pt.openFile(self.hdfActSeqFile,mode="w")
        hdfMakeActiveSeqInd(db=self.dbSql,hdfFile=hdfFile,hdfPath=self.hdfActSeqInd)
        hdfFile.close()

    def reportStat(self):
        """Output a report that shows a high-level overview of data."""
        db.executeAndPrint("""
        select level,avg(cnt) from
        (select level,partaxid,count(*) as cnt from taxa_level group by level,partaxid) a
        group by level
        """)
        db.executeAndPrint("""
        select level,count(*) from
        (select level,partaxid,count(*) as cnt from taxa_level group by level,partaxid) a
        group by level
        """)
        
        ## Some tables to report statistics on data

        ## Group into a single taxid entry per source database
        
        db.createTableAs("src_db_stat","""
        (SELECT src_db        ,
                        taxid ,
                        COUNT(*) AS cnt
                FROM    all_src
                GROUP BY src_db,
                        taxid
        )
        """)

        db.executeAndPrint("""
        SELECT  src_db,
                COUNT(*)
        FROM    src_db_stat
        GROUP BY src_db
        """)

        db.createTableAs("taxa_stat","""
        (SELECT src_type     ,
                        stage,
                        taxid,
                        cat  ,
                        COUNT(*) AS cnt
                FROM    all_src
                GROUP BY src_type,
                        stage    ,
                        taxid    ,
                        cat
        )
        """)
        
        db.executeAndPrint("""
        SELECT  stage,
                COUNT(*)
        FROM    taxa_stat
        GROUP BY stage
        """)

        db.executeAndPrint("""
        SELECT  stage,
                cat  ,
                COUNT(*)
        FROM    taxa_stat
        GROUP BY stage,
                cat
        """)
        
        db.executeAndPrint("""
        SELECT  src_type,
                cat     ,
                COUNT(*)
        FROM    taxa_stat
        GROUP BY src_type,
                cat
        """)
        
        db.executeAndPrint("""
        SELECT  COUNT(*)
        FROM    all_src
        WHERE   kind = 'NC'
        """)

        ## Human must be in NC_ only
        db.executeAndPrint("""
        SELECT  *
        FROM    all_src
        WHERE   taxid = 9606
        """)
        
        db.executeAndPrint("""
        SELECT  rank,
                COUNT(*)
        FROM    all_src
        GROUP BY rank
        """)
        
        db.executeAndPrint("""
        SELECT  stage,
                rank,
                COUNT(*)
        FROM    all_src
        GROUP BY stage, rank
        """)

        db.executeAndPrint("""
        SELECT  divid,
                COUNT(*)
        FROM    seq_sel
        GROUP BY divid
        """)

        db.executeAndPrint("""
        SELECT  rank,
                COUNT(*)
        FROM    seq_sel
        GROUP BY rank
        """)


refseqAccFormat = \
(
('','','','','Undefined'),
('AC','AC_123456','Genomic','Mixed','Alternate complete genomic molecule. This prefix is used for records that are provided to reflect an alternate assembly or annotation. Primarily used for viral, prokaryotic records.'),
('AP','AP_123456','Protein','Mixed','Protein products; alternate protein record. This prefix is used for records that are provided to reflect an alternate assembly or annotation. The AP_ prefix was originally designated for bacterial proteins but this usage was changed.'),
('NC','NC_123456','Genomic','Mixed','Complete genomic molecules including genomes, chromosomes, organelles, plasmids.'),
('NG','NG_123456','Genomic','Mixed','Incomplete genomic region; supplied to support the NCBI genome annotation pipeline. Represents either non-transcribed pseudogenes, or larger regions representing a gene cluster that is difficult to annotate via automatic methods.'),
('NM','NM_123456 NM_123456789','mRNA','Mixed','Transcript products; mature messenger RNA (mRNA) transcripts.'),
('NP','NP_123456 NP_123456789','Protein','Mixed','Protein products; primarily full-length precursor products but may include some partial proteins and mature peptide products.'),
('NR','NR_123456','RNA','Mixed','Non-coding transcripts including structural RNAs, transcribed pseudogenes, and others.'),
('NT','NT_123456','Genomic','Automated','Intermediate genomic assemblies of BAC and/or Whole Genome Shotgun sequence data.'),
('NW','NW_123456 NW_123456789','Genomic','Automated','Intermediate genomic assemblies of BAC or Whole Genome Shotgun sequence data.'),
('NZ','NZ_ABCD12345678','Genomic','Automated','A collection of whole genome shotgun sequence data for a project. Accessions are not tracked between releases. The first four characters following the underscore (e.g. ''ABCD'') identifies a genome project.'),
('XM','XM_123456 XM_123456789','mRNA','Automated','Transcript products; model mRNA provided by a genome annotation process; sequence corresponds to the genomic contig.'),
('XP','XP_123456 XP_123456789','Protein','Automated','Protein products; model proteins provided by a genome annotation process; sequence corresponds to the genomic contig.'),
('XR','XR_123456','RNA','Automated','Transcript products; model non-coding transcripts provided by a genome annotation process; sequence corresponds to the genomic contig.'),
('YP','YP_123456 YP_123456789','Protein','Mixed','Protein products; no corresponding transcript record provided. Primarily used for bacterial, viral, and mitochondrial records.'),
('ZP','ZP_12345678','Protein','Automated','Protein products; annotated on NZ_ accessions (often via computational methods).')
)

class OldTaxaCollector:
    def mergeSelWithSeq(self,skipSeq=False):
        #from itertool import izip
        outFasta = gzip.open(self.selFastaFile,'w',compresslevel=4)
        inpDump = open(self.selDumpFile,'r')
        selGiFile = os.path.abspath(self.selGiFile)
        fldsDump = "gi,taxid,src_db,kind,project,cat,stage,src_type,id,seq_len,divid,rank".split(',')
        pipe = Popen(("fastacmd -i %s -d %s" % (selGiFile,self.srcDbNameAlias)).split(), 
                     cwd=self.blastDataDir, env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        #inpFasta = readFastaRecords(pipe,readSeq=True)
        FGI_SKIP  = 0x01
        FGI_WRITE = 0x02
        FGI_MISM  = 0x04
        giSeen = {}
        iRec = 0
        skip = True
        mismatchRun = 0 # how many gi mismatches in a row
        for line in pipe:
            try:
                if line.startswith(">"):
                    skip = False
                    #header line can be:
                    #>gi|23455713|ref|NC_004301.1| Enterobacteria phage FI, complete genome >gi|15183|emb|X07489.1| Bacteriophage SP genomic RNA
                    headers = line.split('>')
                    gi2 = -1
                    if len(headers) > 2:
                        hdr2 = headers[2]
                        if hdr2.startswith('gi|'):
                            gi2 = int(hdr2.split('|',2)[1])
                            if not giSeen.has_key(gi2):
                                giSeen[gi2] = 0
                    (gifld,gi,accfld,acc,txt) = line[1:].split('|',4)
                    assert gifld == 'gi'
                    gi = int(gi)
                    valsDump = inpDump.readline().rstrip('\n').split(' ') #empty string fields will be ok
                    giDump = int(valsDump[0])
                    taxidDump = int(valsDump[1])
                    try:
                        lineage = self.taxaTree.getNode(taxidDump).lineageRanksStr()
                    except KeyError:
                        lineage = 'NULL'
                        print "Warning: Lineage not found for taxid %s" % (taxidDump,)
                    if giDump == gi:
                        line = ">gi|%s|%s|%s|" % (gi,accfld,acc) + \
                            ''.join(["%s:%s " % (fld,val) for (fld,val) in zip(fldsDump[1:],valsDump[1:])]) + \
                            "lineage:%s " % (lineage,) + \
                            txt
                        if gi2 > 0:
                            giSeen[gi2] |= FGI_WRITE
                        mismatchRun = 0
                    elif giDump == gi2:
                        giSeen[giDump] |= FGI_SKIP
                        skip = True
                        mismatchRun = 0
                    else:
                        if mismatchRun >= 10:
                            raise ValueError("Mismatch between FASTA and SQl DUMP input streams:" + \
                                "fastaTitle = %s valsDump = %s" % (line,' '.join(valsDump)))
                        else:
                            if not giSeen.has_key(giDump):
                                giSeen[giDump] = 0
                            giSeen[giDump] |= (FGI_SKIP | FGI_MISM)
                            skip = True
                            mismatchRun += 1
                            print "GI mismatch for ", giDump
                    if iRec % 10000 == 0:
                        print "Done %s records" % (iRec,)
                    iRec += 1
                elif skipSeq:
                    skip = True
                if not skip:
                    outFasta.write(line)
            except:
                print "Exception with input line: ", line
                pipe.close()
                raise
        pipe.close()
        inpDump.close()
        outFasta.close()
        print 'giSeen = \n', sorted(giSeen.items())
        print 'len(giSeen) = ',len(giSeen)
        for (gi,val) in sorted(giSeen.items()):
            if val & FGI_SKIP and not val & FGI_WRITE:
                print "%s never written %s" % (gi,val)
            if val & FGI_MISM:
                print "%s mistamtch %s" % (gi,val)
    
