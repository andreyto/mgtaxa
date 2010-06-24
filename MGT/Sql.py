### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Taxa import *
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
            ##- it always returns the same value
            self.start = time.time()
            print dedent(sql)

    def __call__(self):
        if self.debug > 0:
            finish = time.time()
            print "SQL finished in %.3f sec" % (finish-self.start)
            self.start = finish

class DbSql(MGTOptions):
    
    def __init__(self):
        MGTOptions.__init__(self)
        atexit.register(dbClose, dbObj=self)
        self.debug = 1

        # This will be lazy constructed because the needed DBAPI module is
        # only available after descendant class __init__ is called
        
        self.numpyTypeMap = None

    def dbapi(self):
        return self.dbmod

    def getNumpyTypeMap(self):
        if self.numpyTypeMap is None:
            self.numpyTypeMap = SqlNumpyTypeMap(self.dbapi())
        return self.numpyTypeMap

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
 
    def dumpCursor(self,cursor):
        pass
    
    def executeAndPrint(self,sql,**kw):
        curs = self.execute(sql,**kw)
        if curs is not None:
            print curs.fetchall()
            curs.close()
    
    def executeAndAssert(self,sql,expect,**kw):
        if self.debug > 0:
            print "Asserting that query will return this: %s" % (expect,)
        curs = self.execute(sql,**kw)
        if curs is None:
            raise SqlResultAssertionError(receive=None,expect=expect)
        rows = curs.fetchall()
        if len(rows) != len(expect):
            raise SqlResultAssertionError(receive=rows,expect=expect)
        for rowReceive, rowExpect in zip(rows,expect):
            if list(rowReceive) != list(rowExpect):
                raise SqlResultAssertionError(receive=rows,expect=expect)
        curs.close()
 
    def executeAndAssertEmpty(self,sql,**kw):
        self.executeAndAssert(sql,tuple(),**kw)
 
    def executeAndAssertZero(self,sql,**kw):
        self.executeAndAssert(sql,((0,),),**kw)

    def selectScalar(self,sql,**kw):
        """Execute sql that must return a single row with a single column and return result as scalar value."""
        curs = self.execute(sql,**kw)
        ret = curs.fetchall()
        assert len(ret) == 1 and len(ret[0]) == 1,"Non-scalar value obtained in 'selectScalar()'"
        ret = ret[0][0]
        curs.close()
        return ret

    def selectAsNx1Dict(self,sql,**kw):
        """Execute sql that must return two columns with Nx1 relation and return result as dict(first->second)."""
        curs = self.execute(sql,**kw)
        ret = curs.fetchall()
        assert len(ret) == 0 or len(ret[0]) == 2,"Result set must be two columns"
        curs.close()
        dret = dict(ret)
        assert len(ret) == len(dret),"Result set must be Nx1 relation. Multi-valued keys were found."
        return dret
    
    def dropTables(self,names):
       for name in names:
           self.dropTable(name)

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

    def createIndices(self,names,table):
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

    def reconnect(self):
        self.close()
        self.open()

    def open(self):
        pass
    
    def close(self):
        pass

    def dialectMatch(self,dialect):
        return dialect is None

    def analyze(self,table):
        self.ddl("ANALYZE TABLE " + table,ifDialect="mysql")

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


    def makeBulkInserterFile(self,**kw):
        k = {}
        k.update(kw)
        k.setdefault("tmpDir",self.tmpDir)
        return BulkInserterFile(self,**k)
   
    def makeBulkInserter(self,*l,**kw):
        return BulkInserter(self,*l,**kw)

    def makeBulkReader(self,*l,**kw):
        return BulkReader(db=self,*l,**kw)

    def saveRecords(self,records,table):
        """Create a new table and save records into it.
        @param records - iterable with each element been a sequence of field values itself
        @param table - instance of SqlTable description class"""
        self.ddl(table.createSql(),dropList=["table "+table.name])
        inserter = self.makeBulkInserterFile(table=table.name)
        for rec in records:
            inserter(rec)
        inserter.flush()

    def createTableFromArray(self,name,arr,withData=True,returnInserter=False,indices=None):
        """Create a table that reflects the fields of Numpy record array.
        @return BulkInserter object or None
        @param name The name of the new table
        @param arr Numpy array to use as template
        @param withData if True, also load data from array
        @param returnInserter if True, return a BulkInserter object; 
        The caller then is resposible for closing the inserter object.
        @param indices, if not None, should be a dictionary with arguments to createIndices, 
        except the 'table' argument, which will be taken from 'name'.
        All fields are constrained as NOT NULL, as Numpy does not have NULL values,
        and NOT NULL constraint speeds up queries.
        """
        numMap = self.getNumpyTypeMap()
        flds = numMap.sqlFromDtype(dtype=arr.dtype)
        table = SqlTable(name=name,fields = [ SqlField(name=f[0],type=f[1],null=False) for f in flds ])
        self.ddl(table.createSql(),dropList=["table "+table.name])
        inserter = None
        if withData or returnInserter:
            inserter = self.makeBulkInserterFile(table=table.name)
        if withData:
            for rec in arr:
                inserter(rec)
            if not returnInserter:
                inserter.close()
        if len(indices) > 0:
            self.createIndices(table=name,**indices)
            self.ddl("analyze table %s" % name,ifDialect="mysql")
        return inserter


    def selectAsArray(self,sql):
        """Execute SQL select and return the entire result set as Numpy record array."""
        reader = self.makeBulkReader(sql=sql)
        ret = reader.allAsArray()
        reader.close()
        return ret

    def exportAsCsv(self,sql,out,withHeader=True):
        """Excecute SQL and export the result as CSV file.
        @todo implement this method by directly saving from result set, w/o the numpy array intermediary."""
        import csv
        reader = self.makeBulkReader(sql=sql,format="list")
        rows = reader.allAsArray()
        names = reader.fieldNames()
        reader.close()
        if isinstance(out,str):
            out = openCompressed(out,'w')
            doClose=True
        else:
            doClose=False
        w = csv.writer(out,delimiter="\t")
        if withHeader:
            w.writerow(names)
        w.writerows(rows)
        if doClose:
            out.close()

def sqlInList(l):
    """Create a string that can be used after SQL 'IN' keyword from a given Python sequence"""
    return "("+','.join(["%s" % x for x in l])+")"


# Classes that describe Sql Db objects as Python objects in a portable way

class SqlField:

    def __init__(self,name,type,null=True):
        self.name = name
        self.type = type
        self.null = null

    def nullSql(self):
        if self.null:
            return ""
        else:
            return "NOT NULL"

class SqlTable:

    def __init__(self,name,fields):
        self.name = name
        self.fields = fields


    def createSql(self):
        return """
        create table %s
        (
        %s
        )
        """ % (self.name,',\n'.join( [ "%s %s %s" % (field.name,field.type,field.nullSql()) for field in self.fields ] ))


class DbSqlLite(DbSql):
    
    def __init__(self,dbpath,strType=str,dryRun=False):
        from pysqlite2 import dbapi2 as dbmod
        self.dbmod = dbmod
        DbSql.__init__(self)
        self.strType = strType
        self.dbpath = dbpath
        self.dryRun = dryRun
        self.con = self.dbmod.connect(self.dbpath)
        self.con.text_factory = self.strType


    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()


class DbSqlMy(DbSql):
    
    def __init__(self,dryRun=False,**kw):
        import MySQLdb as dbmod
        self.dbmod = dbmod
        DbSql.__init__(self)
        self.open(**kw)
        #self.dbmod.server_init(("phyla","--defaults-file=my.cnf"),('server',))


        ## TODO: handle situation due to long periods of computations w/o SQL calls:
        ## Exception _mysql_exceptions.OperationalError: (2006, 'MySQL server has gone away')
        ## or
        ## _mysql_exceptions.OperationalError: (2013, 'Lost connection to MySQL server during query')

    def open(self,**kw):
        import MySQLdb.cursors
        self.close()
        self.con = self.dbmod.connect(port=13306,
                                      #unix_socket="/tmp/atovtchi.mysql.sock",
                                      host=kw.get("host","localhost"),
                                      db=kw.get('db',"mgtaxa"),
                                      user="root",
                                      passwd="OrangeN0",
                                      read_default_file="my.cnf",
                                      read_default_group="client")
                                      ## SSCursor fetches results row-by-row from the server,
                                      ## although I did not see any difference in overall speed
                                      #cursorclass=MySQLdb.cursors.SSCursor)

    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()
            del self.con
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

    def createIndices(self,table,names=None,primary=None,compounds=None,attrib={}):
        """This is a specialization for MySQL, which supports ALTER TABLE ... ADD INDEX ... ADD INDEX.
        @param attrib - optional index attributes. Currently supported is 'unique', e.g. attrib={'id':{'unique':True}}"""
        dropNames = []
        sql = "ALTER TABLE %s " % (table,)
        comma = ""
        if primary is not None:
            sql = sql + "%s\nADD PRIMARY KEY (%s)" % (comma,primary)
            comma = ","
        if names is not None:
            addindex = []
            for name in names:
                dropNames.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                        pass
                addindex.append("ADD %s INDEX ix_%s_%s (%s)" % (unique,table,name,name))
            if len(addindex) > 0:
                sql = sql + "%s\n" % (comma,)+ ",\n".join(addindex)
                comma = ","
        if compounds is not None:
            addindex = []
            for name in compounds.keys():
                dropNames.append(name)
                unique = ""
                try:
                    if attrib[name]["unique"]:
                        unique = "UNIQUE"
                except KeyError:
                        pass
                addindex.append("ADD %s INDEX ix_%s_%s (%s)" % (unique,table,name,compounds[name]))
            if len(addindex) > 0:
                sql = sql + "%s\n" % (comma,) + ",\n".join(addindex)
                comma = ","
        for name in dropNames:
            self.dropIndex(name,table)
        self.ddl(sql)

    def exportToFile(self,sql1,sql2,fileName,fieldsTerm='|',linesTerm=r'\n'):
        rmf(fileName)
        self.ddl("""
        %s
        INTO    OUTFILE '%s'
        FIELDS TERMINATED BY '%s'
        LINES TERMINATED BY '%s'
        %s
        """ % (sql1,fileName,fieldsTerm,linesTerm,sql2))
        

    def exportToStream(self,sql1,sql2,fieldsTerm='|',linesTerm=r'\n'):
        (fobj,bulkFile) = makeTmpFile(dir=self.tmpDir,prefix="sqlBulkExp_",suffix=".csv",createParents=True)
        fobj.close()
        self.exportToFile(sql1=sql1,sql2=sql2,fileName=bulkFile,fieldsTerm=fieldsTerm,linesTerm=linesTerm)
        fobj = open(bulkFile,'r',2**20)
        os.remove(bulkFile) # will cause exception on Windows.
        return fobj


class DbSqlMonet(DbSql):
    
#sql>CREATE USER "root" WITH PASSWORD 'OrangeN0' NAME 'Main User' SCHEMA "sys";
#sql>CREATE SCHEMA "mgtaxa" AUTHORIZATION "root";
#sql>ALTER USER "root" SET SCHEMA "mgtaxa";
	
    def __init__(self,dryRun=False):
        import MonetSQLdb as dbmod
        self.dbmod = dbmod
        DbSql.__init__(self)
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

    def createIndices(self,names,table):
        pass

def createDbSql():
    #db = DbSqlLite("/export/atovtchi/test_seq.db")
    #db = DbSqlMonet()
    db = DbSqlMy()
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


class BulkInserterFile(object):
    
    def __init__(self,db,table,bufLen=500000,workDir="."):
        """Bulk loading from files is not part of SQL standard.
        However, every DBMS seems to have a way to do it, and it is very fast
        (100x for MonetDB) compared to plain 'insert' from 'executemany'.
        This implementation should work for MonetDB and possibly PostgreSQL"""
        from cStringIO import StringIO
        self.db = db
        self.bufLen = bufLen
        self.table = table
        self.workDir = workDir
        self.buf = StringIO()
        self.n = 0
        self.bulkFile = None
        self.isClosed = False

            
    def __call__(self,record):
        buf = self.buf
        escaped = []
        for x in record:
            if isinstance(x,basestring):
                escaped.append('"%s"' % (x,))
            elif x is None: #NULL
                escaped.append('\N') #at least this is how MySQL encodes NULL
            # %s would print bool as True or False, which MySQL will not understand
            elif isinstance(x,(bool,n.bool_)):
                escaped.append("%i" % x)
            else:
                escaped.append(str(x))
        self.buf.write('|'.join(escaped)+'\n')
        self.n += 1
        if self.n == self.bufLen:
            self.flush()

    def __del__(self):
        self.close()

    def newBulkFile(self):
        self.delBulkFile()
        (fobj,self.bulkFile) = makeTmpFile(dir=self.workDir,
                prefix="sqlBulkIns_"+self.table,
                suffix=".csv",
                createParents=True,
                bufsize=0)
        return fobj

    def delBulkFile(self):
        if self.bulkFile is not None:
            os.remove(self.bulkFile)
            self.bulkFile = None

    def close(self):
        if not self.isClosed:
            self.flush()
            self.db = None
            self.delBulkFile()
            self.isClosed = True

    def dumpToBulkFile(self,s):
        # we create a new tmp file every time to circumvent NFS client caching
        # issues that are possible when SQL server and this process are on different hosts.
        # with a new unique file name, worst case is that the server will not find the file,
        # rather than see the old data again.
        out = self.newBulkFile()
        out.write(s)
        fileSync(out)
        out.close()
            
    def flush(self):
        if self.n > 0:
            #self.buf.seek(0)
            #if not hasattr(self,'prev_accum'):
            #    self.prev_accum = debugTaxaCount(self.db)

            s = self.buf.getvalue()
            
            self.dumpToBulkFile(s)

            #out = open(self.bulkFile+'.append','a')
            #out.write(s)
            #out.close()

            self.buf.close()
            self.buf = StringIO()

            self.db.ddl("copy into %s from '%s'" % (self.table,self.bulkFile) )

            #accum = debugTaxaCount(self.db)
            #print "Flushed BulkInserter, now we have %s taxids" % (accum,)
            #assert accum >= self.prev_accum
            #self.prev_accum = accum
            #pdb.set_trace()

            self.n = 0


def debugTaxaCount(db):
    cursor = db.execute("""
    select count(*) from (select taxid from seq group by taxid) a
    """)
    res = cursor.fetchall()
    cursor.close()
    return int(res[0][0])
    

class BulkInserterFileMy(BulkInserterFile):
    """@todo MySQL either has to use 'load data local', with 'local' being enabled both at client and server config,
    or file must be world-readable server-side, and accessible from the server. See MySQL reference for the security
    implications of these alternatives."""

    def __init__(self,*l,**kw):
        BulkInserterFile.__init__(self,*l,**kw)
    
    def newBulkFile(self):
        fobj = super(BulkInserterFileMy,self).newBulkFile()
        #hopefully, it is safe to use chmod on an open file
        chmod(self.bulkFile,'o+r')
        return fobj

    def flush(self):
        if self.n > 0:
            outStr = self.buf.getvalue()
            self.dumpToBulkFile(outStr)
            sql = """
                    load data infile '%s'
                    into table %s
                    fields
                    terminated by '|'
                    optionally enclosed by '"'
                """ % (self.bulkFile,self.table)
            self.db.ddl(sql)
            self.n = 0
            self.buf.close()
            self.buf = StringIO()


class SqlNumpyTypeMap:
    """Map between SQL data types and NumPy data types.
    Precision in digits vs number of bytes for integer NUMBER types are taken
    from MySQL reference (but should be universal)
    http://dev.mysql.com/doc/refman/5.0/en/numeric-types.html
    and from SQL code:
    db.ddl("create table tmp_types (f_bool bool, f_tinyint tinyint, f_smallint smallint, f_med mediumint, f_int int,"+
    " f_big bigint, f_char char(5), f_float float)",dropList=["table tmp_types"])
    curs = db.execute("select * from tmp_types limit 1")
    print curs.description
    """
    
    intDigitsBytes = numpy.array([('bool',1,1,1), ('tinyint',4,1,1), ('smallint',6,2,2),
                                  ('mediumint',9,3,4), ('int',11,4,4), ('bigint',20,8,8)],
                    dtype=[('typeSql', 'S15'), ('digits', 'i4'), ('bytesSql', 'i4'), ('bytesNpy', 'i4')])

    def __init__(self,dbapi):
        """@param dbapi - Python DBAPI module instance."""
        self.dbapi = dbapi

    def close(self):
        self.dbapi = None

    def digitsToBytes(self,dig):
        return self.intDigitsBytes['bytesNpy'][numpy.digitize([dig-0.0001],self.intDigitsBytes['digits'])[0]]
    
    def dtype(self,descr):
        """Take cursor.description object and return NumPy dtype suitable for record array construction.
        dtype is returned in a 'list of tuples' representation e.g. [('id','int8'),('taxid','int4')]
        SQL is case insensitive, so we convert field names to lower case when constructing dtype.
        @bug 'smallint' in MySQL DBAPI has typecode 2, which is not recognized as NUMBER by that DBAPI.
        You should cast 'smallint' fields to some other integer type in your SQL SELECT statement."""

        ## mnemonic names for indexes into cursor.description list
        
        I_NAME  = 0
        I_TYPE  = 1
        I_SIZE  = 3
        I_PREC  = 4
        I_SCALE = 5

        dbapi = self.dbapi
        dt = []
        for fld in descr:
            npy_t = None
            fld_t = fld[I_TYPE]
            if fld_t == dbapi.NUMBER:
                if fld[I_SCALE] > 0:
                    npy_t = 'f8'
                else:
                    npy_t = 'i%s' % (self.digitsToBytes(fld[I_PREC]),)
            elif fld_t == dbapi.STRING:
                npy_t = 'S%s' % (fld[I_SIZE],)
            elif fld_t == dbapi.DATE or fld_t == dbapi.DATETIME:
                ## DATE is string for numpy
                npy_t = 'S%s' % (fld[I_SIZE],)
            else:
                # just try 'double' a a last resort. MySQL dbapi lacks NUMERIC (DECIMAL) constant
                npy_t = 'f8'
                #raise TypeError("Unsuported SQL DBAPI typecode: %s. Field is %s" % (fld_t,fld))
            dt.append((fld[I_NAME].lower(),npy_t))
        return dt

    def sqlFromDtype(self,dtype):
        sql = []
        for fName in dtype.names:
            fDtype = dtype.fields[fName][0]
            kind = fDtype.kind
            itemsize = fDtype.itemsize
            if kind == 'b':
                spec = 'bool'
            elif kind in 'iu':
                if itemsize <= 4:
                    spec = 'integer'
                else:
                    spec = 'bigint'
            elif kind == 'f':
                spec = 'float(%s)' % itemsize
            elif kind in 'cS':
                spec = 'char(%s)' % itemsize
            else:
                raise ValueError("Unsupported Numpy datatype for SQL conversion: %s" % fDtype)
            sql.append((fName,spec))
            #pdb.set_trace()
        return sql


class BulkReader:
    """Class that executes SQL statement and provides an iterator to read results back in chunks as NumPy record arrays or nested lists."""

    def __init__(self,db,sql,bufLen=None,format="array"):
        self.db = db
        self.sql = sql
        curs = db.execute(sql=sql)
        self.curs = curs
        if bufLen is not None:
            curs.arraysize = bufLen
        else:
            bufLen = curs.arraysize
        self.bufLen = bufLen
        self.format = format
        if self.format == "array":
            self.dt = db.getNumpyTypeMap().dtype(curs.description)
        else:
            self.dt = None
        #print "descr = ", curs.description
        #print "dt = ", self.dt
        self.nrows_fetched = 0

    def nrows(self):
        """Return total number of rows in the result, or None if it cannot be determined.
        It looks like not every DBAPI module (or SQL backend) is capable of returning a rowcount."""
        rowcount = self.curs.rowcount
        if rowcount is None or rowcount < 0:
            return None
        else:
            return rowcount

    def nrowsFetched(self):
        return self.nrows_fetched

    def chunks(self):
        """Iterate through result rows in chunks (as numpy arrays or nested lists). 
        Each chunk is of length limited by 'bufLen' argument supplied to the __init__()."""
        curs = self.curs
        dt = self.dt
        format = self.format
        while True:
            rows = curs.fetchmany()
            if not rows:
                break
            #print rows
            self.nrows_fetched += len(rows)
            if format == "array":
                yield numpy.rec.fromrecords(list(rows),dtype=dt)
            else:
                yield rows

    def allAsArray(self):
        """Return all (remaining) records as one numpy record array or nested list."""
        rows = self.curs.fetchall()
        if self.format == "array":
            return numpy.rec.fromrecords(list(rows),dtype=self.dt)
        else:
            return rows

    def fieldNames(self):
        return [ fld[0] for fld in self.curs.description ]

    def close(self):
        self.curs.close()
        self.db = None


