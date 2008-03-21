from MGT.Common import *
import time

def dbClose(dbObj):
    print "Running 'atexit()' handler"
    dbObj.close()

class SqlResultAssertionError(StandardError):
    
    def __init__(self,receive,expect):
        self.receive = receive
        self.expect = expect
    
    def __str__(self):
        return "SqlResultAssertionError\nReceived:\n%s\n\nExpected:\n%s\n" % (self.receive,self.expect)

class SqlWatch:

    def __init__(self,sql,debug):
        self.debug = debug
        if self.debug > 0:
            ##time.clock() seems to be broken on SuSe 10 x86_64 Python 2.4
            ##- always return the same value
            self.start = time.time()
            print sql

    def __call__(self):
        if self.debug > 0:
            finish = time.time()
            print "SQL finished in %.3f sec" % (finish-self.start)
            self.start = finish

class DbSQL:
    
    def __init__(self):
        atexit.register(dbClose, dbObj=self)
        self.debug = 1

    def ddlIgnoreErr(self,*l,**kw):
        curs = self.cursor()
        try:
            curs.execute(*l,**kw)
        except StandardError, msg:
            print msg
            pass
        curs.close()
        
    def ddl(self,sql,ifDialect=None,dropList=tuple(),ignoreError=False,**kw):
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            for dbobj in dropList:
                words = [ x.strip().lower() for x in dbobj.strip().split() ]
                if words[0] == 'index':
                    assert words[2] == 'on'
                    self.dropIndex(words[1],words[3],rawName=True)
                else:
                    try:
                        dropSql = "drop %s" % (dbobj,)
                        if self.debug > 0:
                            print dropSql
                        curs.execute(dropSql)
                    except StandardError, msg:
                        pass
            watch = SqlWatch(sql,self.debug)
            try:
                curs.execute(sql,**kw)
            except StandardError, msg:
                if ignoreError:
                    print msg
                else:
                    raise
            watch()
            curs.close()
        
    def execute(self,sql,ifDialect=None,**kw):
        #if not sql.strip().startswith("drop"):
        #    return
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            watch = SqlWatch(sql,self.debug)
            curs.execute(sql,**kw)
            watch()
            return curs
        else:
            return None

    def executemany(self,sql,data,ifDialect=None,**kw):
        if self.dialectMatch(ifDialect):
            curs = self.cursor()
            watch = SqlWatch(sql,self.debug)
            curs.executemany(sql,data,**kw)
            watch()
            curs.close()
 
    def executeAndPrint(self,sql,**kw):
        curs = self.execute(sql,**kw)
        if curs is not None:
            print curs.fetchall()
            curs.close()
    
    def executeAndAssert(self,sql,expect,**kw):
        curs = self.execute(sql,**kw)
        if curs is not None:
            rows = curs.fetchall()
            if len(rows) != len(expect):
                raise SqlResultAssertionError(receive=rows,expect=expect)
                for rowReceive, rowExpect in zip(rows,expect):
                    if list(rowReceive) != list(rowExpect):
                        raise SqlResultAssertionError(receive=rows,expect=expect)
            curs.close()
 
 
    def dropTable(self,name):
        try:
            self.ddl("drop table " + name)
        except StandardError:
            pass

    def dropIndex(self,name,table,rawName=False):
        """We always create and drop indices as tablename_indexname because in some DBMS (MonetDB) index names should be globally unique."""
        try:
            if rawName:
                self.ddl("drop index %s on %s" % (name,table))
            else:
                self.ddl("drop index ix_%s_%s on %s" % (table,name,table))
        except StandardError:
            pass

    def createColumnIndices(self,names,table):
        """We always create and drop indices as tablename_indexname because in some DBMS (MonetDB) index names should be globally unique.
        This can be specialized for MySQL and others that support ALTER TABLE ... ADD INDEX ... ADD INDEX"""
        for name in names:
            self.dropIndex(name,table)
            self.ddl("create index ix_%s_%s on %s ( %s )" % (table,name,table,name))

    def connection(self):
        return self.con
    
    def cursor(self):
        return self.con.cursor()
    
    def commit(self):
        if hasattr(self,'con'):
            self.con.commit()

    def close(self):
        pass

    def dialectMatch(self,dialect):
        return dialect is None

    def createTableAs(self,name,select):
        """Insulate 'create table as (select ...) [order by ...] with data' from SQL dialect differences.
        @param name - name of table to (re-)create
        @param select - everything that should go between 'create table as' and 'with data'
        Rational: MonetDB Feb2008 requires 'with data' at the end, MySQL 5 does not recognize 'with data'"""

        self.ddl("""
        CREATE TABLE %s AS
        """ % (name,) + select + """
        WITH DATA
        """,
        dropList=["table %s" % (name,)])
        if self.debug:
            curs = self.execute("select count(*) from %s" % (name,))
            if curs is not None:
                print "%s rows created in table %s" % (curs.fetchall(),name)
                curs.close()


    def makeBulkInserterFile(self,*l,**kw):
        return BulkInserterFile(self,*l,**kw)
   
    def makeBulkInserter(self,*l,**kw):
        return BulkInserter(self,*l,**kw)


class DbSQLLite(DbSQL):
    
    def __init__(self,dbpath,strType=str,dryRun=False):
        from pysqlite2 import dbapi2 as dbmod
        self.dbmod = dbmod
        DbSQL.__init__(self)
        self.strType = strType
        self.dbpath = dbpath
        self.dryRun = dryRun
        self.con = self.dbmod.connect(self.dbpath)
        self.con.text_factory = self.strType


    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()


class DbSQLMy(DbSQL):
    
    def __init__(self,dryRun=False):
        import MySQLdb as dbmod
        self.dbmod = dbmod
        DbSQL.__init__(self)        
        #self.dbmod.server_init(("phyla","--defaults-file=my.cnf"),('server',))
        self.con = self.dbmod.connect(unix_socket="/tmp/atovtchi.mysql.sock",
                                      host="localhost",
                                      db="mgtaxa",
                                      user="root",
                                      passwd="OrangeN0",
                                      read_default_file="my.cnf",
                                      read_default_group="client")


    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()
        #self.dbmod.server_end()

    def dropTable(self,name):
        self.ddl("drop table if exists " + name)

    def dialectMatch(self,dialect):
        return dialect is None or dialect == "mysql"

    def createTableAs(self,name,select):
        """Specialization.
        MySQL does not recognize 'with data' suffix"""

        self.ddl("""
        CREATE TABLE %s AS
        """ % (name,) + select + """
        """,
        dropList=["table %s" % (name,)])
        if self.debug:
            curs = self.execute("select count(*) from %s" % (name,))
            if curs is not None:
                print "%s rows created in table %s" % (curs.fetchall(),name)
                curs.close()

    def makeBulkInserterFile(self,*l,**kw):
        return BulkInserterFileMy(self,*l,**kw)


class DbSQLMonet(DbSQL):
    
#sql>CREATE USER "root" WITH PASSWORD 'OrangeN0' NAME 'Main User' SCHEMA "sys";
#sql>CREATE SCHEMA "mgtaxa" AUTHORIZATION "root";
#sql>ALTER USER "root" SET SCHEMA "mgtaxa";
	
    def __init__(self,dryRun=False):
        import MonetSQLdb as dbmod
        self.dbmod = dbmod
        DbSQL.__init__(self)        
        #self.dbmod.server_init(("phyla","--defaults-file=my.cnf"),('server',))
        self.con = self.dbmod.connect(host = 'localhost',
            dbname = 'mgtaxa')#,
            #user = 'root',
            #password = 'OrangeN0')
        return
        self.ddl("""
        CREATE USER 'mgtuser'
        WITH PASSWORD 'OrangeN0'
        NAME 'Power user of MGTAXA'
        SCHEMA 'sys'
        """,
        ignoreError=True)
        self.ddl("""
        CREATE SCHEMA "mgtschema" AUTHORIZATION "mgtuser";
        """,
        ignoreError=True)
        self.ddl("""
        ALTER USER "mgtuser" SET SCHEMA "mgtschema"
        """,
        ignoreError=True)
        self.con.close()
        self.con = self.dbmod.connect(host = 'localhost',
            dbname = 'mgtaxa',
            lang = 'sql',
            user = 'mgtuser',
            password = 'OrangeN0')
        

    def close(self):
        if hasattr(self,'con'):
            self.con.close()
        #self.dbmod.server_end()

    def dropIndex(self,name,table,rawName=False):
        """MonetDB does not recognize the 'on <table name>' clause in 'drop index' command'"""
        try:
            if rawName:
                self.ddl("drop index %s" % (name,))
            else:
                self.ddl("drop index ix_%s_%s" % (table,name))
        except StandardError:
            pass

    def createColumnIndices(self,names,table):
        pass

def createDbSQL():
    #db = DbSQLLite("/export/atovtchi/test_seq.db")
    db = DbSQLMonet()
    #db = DbSQLMy()
    return db


class IntIdGenerator(object):
    def __init__(self,start=1):
        self.n = start - 1

    def __call__(self):
        self.n += 1
        return self.n

    def last(self):
        return self.n


class BulkInserter:
    
    def __init__(self,db,sql,bufLen):
        self.db = db
        self.bufLen = bufLen
        self.sql = sql
        self.buf = []
        
    def __call__(self,record):
        self.buf.append(record)
        if len(self.buf) == self.bufLen:
            self.flush()

    def __del__(self):
        self.db = None
            
    def flush(self):
        if len(self.buf) > 0:
            self.db.executemany(self.sql,self.buf)
            self.buf = []


class BulkInserterFile:
    
    def __init__(self,db,sql,bufLen):
        """Bulk loading from files is not part of SQL standard.
        However, every DBMS seems to have a way to do it, and it is very fast
        (100x for MonetDB) compared to plain 'insert' from 'executemany'.
        This implementation should work for MonetDB and possibly PostgreSQL"""
        from cStringIO import StringIO
        self.db = db
        self.bufLen = bufLen
        self.sql = sql
        self.buf = StringIO()
        self.n = 0
        self.bulkFile = "/usr/local/scratch/atovtchi/bulk.tmp"

    def __call__new(self,record):
        buf = self.buf
        #apparently it's ok for MonetDB to enclose any type in ""
        escaped = '|'.join(['"%s"' % (x,) for x in record])
        self.buf.write(escaped+'\n')
        self.n += 1
        if self.n == self.bufLen:
            self.flush()
            
    def __call__(self,record):
        buf = self.buf
        escaped = []
        for x in record:
            if isinstance(x,basestring):
                escaped.append('"%s"' % (x,))
            else:
                escaped.append(str(x))
        self.buf.write('|'.join(escaped)+'\n')
        self.n += 1
        ##TODO: until we convert this to using staging
        ##table, just call 'flush()' once for the entire
        ##input set
        #if self.n == self.bufLen:
        #    self.flush()

    def __del__(self):
        self.db = None
            
    def flush(self):
        if self.n > 0:
            #self.buf.seek(0)
            #if not hasattr(self,'prev_accum'):
            #    self.prev_accum = debugTaxaCount(self.db)

            s = self.buf.getvalue()
            
            out = open(self.bulkFile,'w')
            out.write(s)
            out.close()

            #out = open(self.bulkFile+'.append','a')
            #out.write(s)
            #out.close()

            self.buf.close()
            self.buf = StringIO()

            self.db.ddl("copy into %s from '%s'" % (self.sql,self.bulkFile) )

            #accum = debugTaxaCount(self.db)
            #print "Flushed BulkInserter, now we have %s taxids" % (accum,)
            #assert accum >= self.prev_accum
            #self.prev_accum = accum
            #pdb.set_trace()
            os.remove(self.bulkFile)

            self.n = 0


def debugTaxaCount(db):
    cursor = db.execute("""
    select count(*) from (select taxid from mgt_seq group by taxid) a
    """)
    res = cursor.fetchall()
    cursor.close()
    return int(res[0][0])
    

class BulkInserterFileMy(BulkInserterFile):
    """@todo MySQL either has to use 'load data local', with 'local' being enabled both at client and server config,
    or file must be world-readable server-side, and accessible from the server. See MySQL reference for the security
    implications of these alternatives."""
    def flush(self):
        if self.n > 0:
            try:
                out = open(self.bulkFile,'w')
                outStr = self.buf.getvalue()
                out.write(outStr)
                print outStr[:200]
                out.close()
                chmod(self.bulkFile,'o+r')
                sql = """
                        load data infile '%s'
                        into table %s
                        fields
                        terminated by '|'
                        optionally enclosed by '"'
                    """ % (self.bulkFile,self.sql)
                self.db.ddl(sql)
                self.n = 0
                self.buf.close()
                self.buf = StringIO()
            finally:
                os.remove(self.bulkFile)


class DbSeqSource(PhyOptions):
    """Database of sequence source for training the classifier"""
    
    def __init__(self,dbSql):
        
        PhyOptions.__init__(self)
        
        self.dbSql = dbSql
        self.ncbiDbs =  (
                         Struct(id='g',db='refseq_genomic'),
                         Struct(id='o',db='other_genomic'),
                         Struct(id='n',db='nt'),
                         Struct(id='h',db='htgs'),
                         Struct(id='w',db='wgs')
                         )


    def makeBlastAlias(self,idList=None,dbNames=None):
        from glob import glob
        from datetime import datetime
        from textwrap import dedent
        if dbNames is None:
            dbNames = [ db.db for db in self.ncbiDbs ]
        directDbRefs = False
        #comment out GILIST field in alias file if idList is None
        if idList is None:
            comGILIST = '#'
        else:
            comGILIST = ''
        dbNameAlias = self.srcDbNameAlias
        idListBin = dbNameAlias+'.gil'
        aliasStr = dedent("""\
        #
        # Alias file created %s
        #
        #
        TITLE Custom Unified DB for Phyla
        #
        DBLIST %%s
        #
        %sGILIST %s
        #
        #OIDLIST
        #
        """ % (datetime.today().ctime(),comGILIST,idListBin))
        cwd = os.getcwd()
        try:
            os.chdir(self.blastDataDir)
            if directDbRefs:
                chunks = []
                for rec in dbNames:
                    chunks += sorted(list(set(
                                    ( y.group(1) for y in  
                                      (re.match('('+rec.db+'\.[0-9]+)\..*',x) 
                                       for x in glob(rec.db+'.*')) 
                                    if y is not None )
                                    )))
            else:
                chunks = dbNames
            strToFile(aliasStr % ' '.join(chunks),dbNameAlias+'.nal')
        finally:
            os.chdir(cwd)
        
    def mergeSelWithSeq(self,skipSeq=False):
        #from itertool import izip
        outFasta = gzip.open(self.selFastaFile,'w',compresslevel=4)
        inpDump = open(self.selDumpFile,'r')
        selGiFile = os.path.abspath(self.selGiFile)
        fldsDump = "gi,taxid,src_db,kind,project,cat,stage,src_type,iid,seq_len,divid,rank".split(',')
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


    def rebuild(self):
        #self.loadGiTaxNumpy()
        self.loadTaxCategories()
        self.loadTaxNodes()
        self.loadRefseqAcc()
        self.loadSeq()
        self.selectTaxSet()

    def loadSeq(self):
        db = self.dbSql
        self.loadGiTaxPickled()
        self.createTableSeq()
        sql = """
                load data infile '/usr/local/scratch/atovtchi.mgt_seq.bulk.tmp'
                into table mgt_seq
                fields
                terminated by '|'
                optionally enclosed by '"'
        """
        sql = """
        copy into mgt_seq from '/usr/local/scratch/atovtchi.mgt_seq.bulk.tmp'
        """
        db.ddl(sql)

        #self.idGenSeq = IntIdGenerator()
        #inserter = db.makeBulkInserterFile(sql='mgt_seq',bufLen=500000)
        #self.loadSeqNCBI(self.ncbiDbs[0],inserter)
        #self.loadSeqNCBI(self.ncbiDbs[1],inserter)
        #self.loadSeqNCBI(self.ncbiDbs[2],inserter)
        #self.loadSeqNCBI(self.ncbiDbs[3],inserter)
        #self.loadSeqNCBI(self.ncbiDbs[4],inserter)
        #inserter.flush()
        db.ddl("analyze table mgt_seq",ifDialect="mysql")
        db.createColumnIndices(names=["taxid","src_db","kind","project"],table="mgt_seq")
        self.clearGiTax()

    def clearGiTax(self):
        self.gi2taxa = None
    
    def createTableSeq(self):
        # we removed "auto_increment primary key" from iid
        # to speed up bulk load
        self.dbSql.ddl("""
        create table mgt_seq
        (
        iid integer,
        gi bigint,
        taxid integer,
        src_db char(1),
        project char(4),
        seq_len bigint,
        acc char(20),
        kind char(2),
        seq_hdr char(%s)
        )
        """ % (self.fastaHdrSqlLen,),
        dropList=["table mgt_seq"])
        
    def loadSeqNCBI(self,db,inserter):
        pipe = Popen(("fastacmd -D 1 -d " + db.db).split(), cwd=self.blastDataDir, env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        #inp = readFastaRecords(pipe,readSeq=False)
        inp = FastaReader(pipe)
        #curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 50000
        sql = """
        insert into mgt_seq
        (iid, gi, taxid, src_db, project, seq_len, acc, kind, seq_hdr)
        values
        (%s,%s,%s,%s,%s,%s,%s,%s,%s)
        """
        #inserter = self.dbSql.makeBulkInserter(sql=sql,bufLen=bufLen)
        #inserter.n = 1
        #inserter.flush()
        #inp.close()
        #return
        for rec in inp.records():
            title = rec.header()[1:-1] #remove '>' and '\n'
            (gifld,gi,accfld,acc,txt) = title.split('|',4)
            assert gifld == 'gi'
            gi = int(gi)
            try:
                taxid = int(self.gi2taxa[gi])
            except KeyError:
                print "Taxid not found for gi: "+gi
                taxid = 0
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
            values = (self.idGenSeq(),gi,taxid,db.id,project,seqLen,acc,kind,title[:self.fastaHdrSqlLen])
            #values = [ str(x) for x in values ]
            inserter(values)
            if iRec % 50000 == 0:
                print db.db, title, taxid, iRec, seqLen
                #if iRec >= 500000:
                #    break
            iRec += 1
        inp.close()

    def loadRefseqAcc(self):
        self.dbSql.ddl("""\
        create table mgt_refseq_acc
        (
        prefix char(2),
        acc char(40),
        molecule char(10),
        method char(15),
        descr varchar(400)
        )
        """,
        dropList=["table mgt_refseq_acc"])
        self.dbSql.createColumnIndices(names=["prefix"],table="mgt_refseq_acc")
        
        self.dbSql.executemany("""\
        insert into mgt_refseq_acc
        (prefix,acc,molecule,method,descr)
        values
        (%s,%s,%s,%s,%s)
        """,refseqAccFormat)
        
    def loadTaxCategories(self):
        self.dbSql.ddl("""\
        create table mgt_taxa_cat
        (cat char(1), taxid_ancestor integer, taxid integer)
        """,
        dropList=["table mgt_taxa_cat"])
        inp = open(self.taxaCatFile,'r')
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 500000
        sql = """\
        insert into mgt_taxa_cat
        (cat,taxid_ancestor,taxid)
        values
        (%s,%s,%s)
        """
        inserter = self.dbSql.makeBulkInserterFile(sql="mgt_taxa_cat",bufLen=bufLen)
        for rec in inp:
            vals = rec.split()
            vals[1] = int(vals[1])
            vals[2] = int(vals[2])
            inserter(vals)
        inp.close()
        inserter.flush()
        curs.close()
        self.dbSql.ddl("""alter table mgt_taxa_cat ADD PRIMARY KEY (taxid)""",ifDialect="mysql")
        self.dbSql.createColumnIndices(names=["cat","taxid_ancestor"],table="mgt_taxa_cat")

    def loadTaxNodes(self):
        self.dbSql.ddl("""\
        create table mgt_taxa_node
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
        dropList=["table mgt_taxa_node"])
        inp = open(self.taxaNodesFile,'r')
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 500000
        sql = """\
        insert into mgt_taxa_node
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
        inserter = self.dbSql.makeBulkInserterFile(sql="mgt_taxa_node",bufLen=bufLen)
        for rec in inp:
            # The last '|' is behind the last field, so we cut it
            inserter([ x.strip() for x in rec.split('|')][:-1])
        inserter.flush()
        inp.close()
        curs.close()
        #Get rid of spaces in 'rank' field to simplify export of data later on
        self.dbSql.execute("update mgt_taxa_node set rank = replace(rank,' ','_')").close()
        self.dbSql.ddl("""alter table mgt_taxa_node ADD PRIMARY KEY (taxid)""",ifDialect="mysql")
        self.dbSql.createColumnIndices(names=["partaxid","divid"],table="mgt_taxa_node")


    def loadTaxNodesMem(self):
        self.taxaTree = TaxaTree(ncbiDumpFile=self.taxaNodesFile)
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
        inserter = self.dbSql.makeBulkInserterFile(sql='gi_taxa',bufLen=bufLen)
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
        inp = openGzip(self.taxaGiFile,'r')
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


    def selectTaxSet(self):
        db = self.dbSql
        db.createTableAs("mgt_taxa_src_1","""
            (SELECT taxid                ,
                            src_db       ,
                            kind         ,
                            project      ,
                            0            AS priority,
                            COUNT(*)     AS cnt     ,
                            SUM(seq_len) AS seq_len
                    FROM    mgt_seq
                    GROUP BY taxid ,
                            src_db ,
                            kind   ,
                            project
            )
        """)

        db.createTableAs("mgt_taxa_src","""
                (SELECT *
                        FROM    mgt_taxa_src_1
                )
        ORDER BY taxid ,
                src_db ,
                kind   ,
                project
        """)
        
        db.dropTable("mgt_taxa_src_1")
        
        #db.ddl("""ALTER TABLE mgt_taxa_src ADD id INTEGER auto_increment PRIMARY KEY""",ifDialect="mysql")
        db.ddl("""ALTER TABLE mgt_taxa_src ADD id INTEGER auto_increment  PRIMARY KEY""")
        db.createColumnIndices(names=["taxid","src_db","kind","project"],table="mgt_taxa_src")
        
        ## Drop "Alternative" full genomes (e.g. haplotypes)
        db.ddl("""
        DELETE
        FROM    mgt_taxa_src
        WHERE   kind = 'AC'
        """)

        db.createTableAs("mgt_wgs_src_1","""
        (SELECT *
                FROM    mgt_taxa_src
                WHERE   project <> ''
                    AND kind    <> ''
        )
        """)
        
        db.ddl("""
        INSERT
        INTO    mgt_wgs_src_1
        SELECT  a.*
        FROM    mgt_taxa_src a
        WHERE   a.project   <> ''
            AND a.taxid NOT IN
                (SELECT taxid
                FROM    mgt_wgs_src_1
                )
        """)

        db.createTableAs("mgt_wgs_src_2","""
        (SELECT a.*
                FROM    mgt_wgs_src_1 a
                WHERE   a.seq_len >=
                        (SELECT MAX(b.seq_len)
                        FROM    mgt_wgs_src_1 b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)

        db.createTableAs("mgt_wgs_src_3","""
        (SELECT a.*
                FROM    mgt_wgs_src_2 a
                WHERE   a.id >=
                        (SELECT MAX(b.id)
                        FROM    mgt_wgs_src_2 b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)
        db.dropTable("mgt_wgs_src_2")

        db.createTableAs("mgt_wgs_src_nr","""
        (SELECT *
                FROM    mgt_wgs_src_3
        )
        """)

        db.createTableAs("mgt_nc_src_1","""
        (SELECT a.*
                FROM    mgt_taxa_src a
                WHERE   a.kind        = 'NC'
                    AND a.seq_len * 2 >
                        (SELECT MAX(b.seq_len)
                        FROM    mgt_wgs_src_nr b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)
        
        db.createTableAs("mgt_htg_over_nc_src_1","""
        (SELECT a.*
                FROM    mgt_taxa_src a,
                        mgt_nc_src_1 b
                WHERE   a.src_db     = 'h'
                    AND a.kind      <> 'NC'
                    AND a.taxid      = b.taxid
                    AND a.seq_len    > b.seq_len * 5
                    AND b.seq_len    < 500000
                    AND a.taxid NOT IN
                        (SELECT taxid
                        FROM    mgt_wgs_src_nr
                        )
        )
        """)

        db.createTableAs("mgt_nc_src_2","""
        (SELECT *
                FROM    mgt_nc_src_1
                WHERE   taxid NOT IN
                        (SELECT taxid
                        FROM    mgt_htg_over_nc_src_1
                        )
        )
        """)
        
        db.createTableAs("mgt_nc_src_3","""
        (SELECT a.*
                FROM    mgt_nc_src_2 a
                WHERE   a.id >=
                        (SELECT MAX(b.id)
                        FROM    mgt_nc_src_2 b
                        WHERE   a.taxid = b.taxid
                        )
        )
        """)
        
        db.createTableAs("mgt_nc_src","""
        (SELECT a.*,
                        'nc' AS stage
                FROM    mgt_nc_src_3 a
        )
        """)
        
        db.createTableAs("mgt_wgs_src","""
        (SELECT a.*,
                        'wg' AS stage
                FROM    mgt_wgs_src_nr a
                WHERE   taxid NOT IN
                        (SELECT taxid
                        FROM    mgt_nc_src
                        )
        )
        """)
        
        db.createTableAs("mgt_htg_src_1","""
        (SELECT *
                FROM    mgt_taxa_src
                WHERE ( src_db      = 'h'
                     OR kind       <> '' )
                    AND (taxid NOT IN
                        (SELECT taxid
                        FROM    mgt_nc_src
                        )
                    AND taxid NOT IN
                        (SELECT taxid
                        FROM    mgt_wgs_src
                        ) )
        )
        """)
        
        db.createTableAs("mgt_htg_src","""
        (SELECT a.*,
                        'ht' AS stage
                FROM    mgt_htg_src_1 a
        )
        """)
        
        db.createTableAs("mgt_gen_src","""
        (SELECT *
                FROM    mgt_nc_src
        )
        """)
        
        db.ddl("""
        INSERT
        INTO    mgt_gen_src
        SELECT  *
        FROM    mgt_wgs_src
        
        UNION
        
        SELECT  *
        FROM    mgt_htg_src
        """)

        db.createColumnIndices(names=["taxid"],table="mgt_gen_src")
        
        db.createTableAs("mgt_nt_src","""
        (SELECT a.*,
                        'nt' AS stage
                FROM    mgt_taxa_src a
                WHERE   taxid NOT IN
                        (SELECT taxid
                        FROM    mgt_gen_src
                        )
        )
        """)

        db.createTableAs("mgt_all_src_1","""
        (SELECT a.*,
                        'g' AS src_type
                FROM    mgt_gen_src a
        )
        """)
        
        db.ddl("""
        INSERT
        INTO    mgt_all_src_1
        SELECT  a.*,
                'o' AS src_type
        FROM    mgt_nt_src a
        """)
    
    #def selectTaxSet(self):

        #db = self.dbSql
        
        db.createTableAs("mgt_all_src","""
        (SELECT a.*            ,
                        b.cat  ,
                        c.divid,
                        c.rank
                FROM    mgt_all_src_1 a
                        LEFT JOIN mgt_taxa_cat b
                        ON      a.taxid = b.taxid
                        LEFT JOIN mgt_taxa_node c
                        ON      a.taxid = c.taxid
        )
        """)

        db.createColumnIndices(names=["taxid","src_db","kind","project","cat","divid","rank"],
            table="mgt_all_src")
        db.ddl("ALTER TABLE mgt_all_src ADD PRIMARY KEY (id)",ifDialect="mysql")
        db.ddl("""CREATE INDEX ix_mgt_all_src_src ON mgt_all_src (taxid,src_db,kind,project)""",
            dropList=["index ix_mgt_all_src_src on mgt_all_src"],ifDialect="mysql")

        ## Some tables to report statistics on data

        ## Group into a single taxid entry per source database
        
        db.createTableAs("mgt_src_db_stat","""
        (SELECT src_db        ,
                        taxid ,
                        COUNT(*) AS cnt
                FROM    mgt_all_src
                GROUP BY src_db,
                        taxid
        )
        """)

        db.executeAndPrint("""
        SELECT  src_db,
                COUNT(*)
        FROM    mgt_src_db_stat
        GROUP BY src_db
        """)

        db.createTableAs("mgt_taxa_stat","""
        (SELECT src_type     ,
                        stage,
                        taxid,
                        cat  ,
                        COUNT(*) AS cnt
                FROM    mgt_all_src
                GROUP BY src_type,
                        stage    ,
                        taxid    ,
                        cat
        )
        """)
        
        db.executeAndPrint("""
        SELECT  stage,
                COUNT(*)
        FROM    mgt_taxa_stat
        GROUP BY stage
        """)

        db.executeAndPrint("""
        SELECT  stage,
                cat  ,
                COUNT(*)
        FROM    mgt_taxa_stat
        GROUP BY stage,
                cat
        """)
        
        db.executeAndPrint("""
        SELECT  src_type,
                cat     ,
                COUNT(*)
        FROM    mgt_taxa_stat
        GROUP BY src_type,
                cat
        """)
        
        db.executeAndPrint("""
        SELECT  COUNT(*)
        FROM    mgt_all_src
        WHERE   kind = 'NC'
        """)

        ## Human must be in NC_ only
        db.executeAndPrint("""
        SELECT  *
        FROM    mgt_all_src
        WHERE   taxid = 9606
        """)

        return

        db.ddl("""CREATE INDEX ix_mgt_seq_src ON mgt_seq (taxid,src_db,kind,project)""",
            dropList=["index ix_mgt_seq_src on seq"],ifDialect="mysql")

        ## in MySQL computation of acc_no_ver (ACC w/o '.X' suffix)
        ## would be just SUBSTRING_INDEX( a.acc , '.', 1 ) AS acc_no_ver,
        ## but we use SQL standard conforming expression
        
        db.createTableAs("mgt_seq_sel","""
        (SELECT a.*                                                    ,
                        b.cat                                          ,
                        b.stage                                        ,
                        b.src_type                                     ,
                        SUBSTRING(a.acc FROM 1 FOR
                        CASE
                           WHEN POSITION('.' IN a.acc) <> 0
                           THEN POSITION('.' IN a.acc) - 1
                           ELSE LENGTH(a.acc)
                        END)
                        AS acc_no_ver                                  ,
                        b.divid                                        ,
                        b.rank
                FROM    mgt_seq a,
                        mgt_all_src b
                WHERE   a.taxid   = b.taxid
                    AND a.src_db  = b.src_db
                    AND a.kind    = b.kind
                    AND a.project = b.project
                    AND a.taxid  <> 0
        )
        """)
        
        db.createColumnIndices(names=["taxid","src_db","kind","project","cat",
            "stage","src_type","gi","acc","acc_no_ver","divid","rank"],
            table="mgt_seq_sel")

        db.ddl("ALTER TABLE mgt_seq_sel ADD PRIMARY KEY (iid)",ifDialect="mysql")

        ## Proved to be essential in MySQL for efficient planning of queries
        db.ddl("analyze table mgt_seq_sel",ifDialect="mysql")

        db.executeAndPrint("""
        SELECT  divid,
                COUNT(*)
        FROM    mgt_seq_sel
        GROUP BY divid
        """)

        db.executeAndPrint("""
        SELECT  rank,
                COUNT(*)
        FROM    mgt_seq_sel
        GROUP BY rank
        """)

        ## Backup all records from seq_sel that are to be discarded in case we
        ## want to analyze them later
        
        db.ddl("""
        CREATE TABLE mgt_seq_sel_del (LIKE seq_sel)
        """,
        dropList=["table mgt_seq_sel_del"])

        db.ddl("ALTER TABLE mgt_seq_sel_del DISABLE KEYS",ifDialect="mysql")
        
        db.ddl("""
        INSERT
        INTO    mgt_seq_sel_del
        SELECT  *
        FROM    mgt_seq_sel
        WHERE   taxid IS NULL
            OR divid IS NULL
            OR cat   IS NULL
            OR divid IN (7,8,11)
        """)
        
        db.ddl("ALTER TABLE mgt_seq_sel_del ENABLE KEYS",ifDialect="mysql")
        
        db.ddl("""
        DELETE
        FROM    mgt_seq_sel
        WHERE   iid IN
                (SELECT iid
                FROM    mgt_seq_sel_del
                )
        """)

        db.ddl("ANALYZE TABLE mgt_seq_sel",ifDialect="mysql")

        return

        ## save gis only into a file to be used as gi list for NCBI alias db
        
        db.ddl("""
        SELECT  gi
        INTO    OUTFILE '/home/atovtchi/scratch/mgtdata/phyla_sel.gi'
        FIELDS TERMINATED BY ' '
        LINES TERMINATED BY '\n'
        FROM    mgt_seq_sel
        ORDER BY taxid,
                iid
        """)

        ## save all importand fields for selected seq records into a
        ## flat file to be merged with
        ## fasta sequence data extracted from blast databases
        ## through just created gi list
        ## the order must match the preceding 'select' for gis
        ## Full csv format will have something like OPTIONALLY ENCLOSED BY '"'

        db.ddl("""
        SELECT  gi      ,
                taxid   ,
                src_db  ,
                kind    ,
                project ,
                cat     ,
                stage   ,
                src_type,
                iid     ,
                seq_len ,
                divid   ,
                rank
        INTO    OUTFILE '/home/atovtchi/scratch/mgtdata/phyla_sel.csv'
        FIELDS TERMINATED BY ' '
        LINES TERMINATED BY '\n'
        FROM    mgt_seq_sel
        ORDER BY taxid,
                iid
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


    
