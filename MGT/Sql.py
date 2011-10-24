### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
import time, csv

def dbClose(dbObj):
    #print "Running 'atexit()' handler"
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

## Classes that describe definitions of SQL fields and table
## These are used primarily to construct DDL statements.
class SqlField:
    """This class defines SQL table field"""

    def __init__(self,name=None,type=None,null=True):
        self.name = name
        self.type = type
        self.null = null

    def nullSql(self):
        if self.null:
            return ""
        else:
            return "NOT NULL"

    def tableDef(self):
        return " ".join( [ (attr if attr else "") for attr in (self.name,self.type, self.nullSql()) ] )

class SqlTable:
    """This class defines SQL table"""

    @classmethod
    def fromDb(klass,db,name):
        """Generate SqlTable instance by querying an actual table in the database.
        @param db DbSql instance
        @param name Table name
        @return new SqlTable instance
        """
        ## mnemonic names for indexes into cursor.description list
        
        I_NAME  = 0
        I_TYPE  = 1
        I_SIZE  = 3
        I_PREC  = 4
        I_SCALE = 5
        
        descr = db.getTableDescr(name)
        fields = [ SqlField(name=f[I_NAME],type="%s(%s)" % (f[I_TYPE],f[I_SIZE])) \
                for f in descr ]
        return klass(name=name,fields=fields)

    def __init__(self,name,fields):
        self.name = name
        self.fields = fields


    def createSql(self):
        """Return SQL DDL string that constructs this table"""
        return """
        create table %s
        (
        %s
        )
        """ % (self.name, ",\n".join( [ field.tableDef() for field in self.fields ] ))

    def insertSql(self):
        """Return SQL DML string that can be passed to Python DB-API cursor.executemany() method"""
        return """
        insert into %s
        (
        %s
        )
        values
        (
        %s
        )
        """ % (self.name, ",\n".join( [ field.name for field in self.fields ]), ",\n".join( [ "?" for field in self.fields ]))


class DbSql(MGTOptions):
    """Wrapper around DB-API Connection class with some convenience methods.
    @todo Convert this to using SQL Alchemy. SQL Alchemy imposed too big abstraction penalty in the past,
    but this might not be the case anymore assuming carefully following its best use practices.
    We will still likely to need our bulk loading methods.
    """
        
    ## Default SQL field definition - a fall-back field type to create
    ## when nothing more specific is provided by the user
    defField = SqlField(name="fld",type="char(40)")
    
    def __init__(self):
        MGTOptions.__init__(self)
        atexit.register(dbClose, dbObj=self)
        self.debug = 0

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

    def getTableDescr(self,name):
        """Return table description as seen by DB API module"""
        curs = self.execute("select * from %s limit 1" % (name,))
        descr = curs.description
        curs.close()
        return descr

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

    def selectAll(self,sql,**kw):
        """Convenience method that does for select statement and execute+fetchall in one step.
        Use for statements with small result set.
        @param sql SQL SELECT statement
        @return result of cursor.fetchall (sequence of tuples)"""
        curs = self.execute(sql,**kw)
        ret = curs.fetchall()
        curs.close()
        return ret
    
    def selectScalar(self,sql,**kw):
        """Execute sql that must return a single row with a single column and return result as scalar value."""
        ret = self.selectAll(sql=sql,**kw)
        assert len(ret) == 1 and len(ret[0]) == 1,"Non-scalar value obtained in 'selectScalar()'"
        ret = ret[0][0]
        return ret

    def selectAsNx1Dict(self,sql,**kw):
        """Execute sql that must return two columns with Nx1 relation and return result as dict(first->second)."""
        ret = self.selectAll(sql=sql,**kw)
        assert len(ret) == 0 or len(ret[0]) == 2,"Result set must be two columns"
        dret = dict(ret)
        assert len(ret) == len(dret),"Result set must be Nx1 relation. Multi-valued keys were found."
        return dret
    
    def selectAs1Col(self,sql,**kw):
        """Execute sql that must return one column and return result as 1D sequence."""
        ret = self.selectAll(sql=sql,**kw)
        assert len(ret) == 0 or len(ret[0]) == 1,"Result set must have one column"
        return [ row[0] for row in ret ]
    
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

    def createTableAs(self,name,select,indices=None):
        """Save the results of SQL SELECT as a new table.
        This abstracts "create ... as ..." operation from minor differences in SQL dialects.
        For example, MonetDB Feb2008 required 'with data' at the end, MySQL 5 and SQLite do 
        not recognize 'with data'
        Override in derived classes if necessary.
        @param name name of table to (re-)create - existing table will be replaced
        @param select SQL select statement
        @param indices if present, will be passed to createIndices() method
        """
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
        if indices is not None and len(indices) > 0:
            self.createIndices(table=name,**indices)
            self.ddl("analyze table %s" % name,ifDialect="mysql")

    def makeBulkInserterFile(self,**kw):
        kw = copy(kw)
        if hasattr(kw["table"],"name"):
            kw["table"] = kw["table"].name
        kw.setdefault("tmpDir",self.tmpDir)
        return BulkInserterFile(self,**kw)
   
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
        if indices is not None and len(indices) > 0:
            self.createIndices(table=name,**indices)
            self.ddl("analyze table %s" % name,ifDialect="mysql")
        return inserter

    def createTableFromCsv(self,name,csvFile,fieldsMap={},defField=None,hasHeader=False,
            dialect="excel-tab",dialect_options={},indices=None,preProc=None,hdrPreProc=None):
        """Create and fill a table from CSV file.
        The intention is to provide a single command to load a CSV file into SQL table
        where reasonable default values can be generated for all options.
        @param name The name of the new table
        @param csv Either a file name, in which case csv.reader(openCompresed(),dialect=dialect) 
        will be used to open the file, or it should be an existing file stream, or csv.reader object
        (in the latter case dialect and dialect_options parameters are ignored).
        @param fieldsMap A dictionary that either maps field position to SqlField instances
        (if hasHeader is False) or maps field names to SqlField instances (if hasHeader is True).
        In the latter case, each SqlField instance can replace the name from the header, or leave the
        name as None, in which case the name supplied by the header will be used.
        @param defField SqlField instance to generate default SQL definitions for fields not
        mapped by fieldsMap. If hasHeader is False, defField.name will be used as a 
        prefix, such that the field name becomes prefix_xxx where xxx is the absolute field 
        position in CSV row. Otherwise, defField.name is ignored, and only the other attributes are
        used to define default field type etc.
        @param hasHeader tells whether to treat the first row of CSV file as header
        @param dialect Dialect string defined in csv module
        @param dialect_options Passed to csv.reader(**dialect_options)
        @param indices, if not None, should be a dictionary with arguments to createIndices, 
        except the 'table' argument, which will be taken from 'name'
        @param preProc If specified, this should be a method that will be applied to each
        row returned by csv.reader. The method must return a sequence (possibly empty) of new
        rows, which will be inserted into SQL table instead of the original row. Their size and
        types must match the original row. The preProc must have this signature:
        preProc(row,fields,nameToInd) where row is returned by csv.reader.next(); fields is a list
        of SqlField objects matching the row fields; nameToInd is a dictionary mapping field names
        to indexes of fields in the row. The last two parameters allow preProc's code to access row
        elements by field names, e.g. row[nameToInd["seqid"]]. The preProc parameter addresses a common
        use case where the input file is very large but we need to load into the SQL DB only a small
        subset of it for which a simple filter condition exists such as set membership. It also covers
        simple manipulation of input data such as various string substitutions.
        Example:
        preProc = lambda row,fields,nameToInd,idSet=set(1,2,3): \
                ( row, ) if row[namesToInd["seqId"]] in idSet else (,)
        createTableFromCsv(...,preProc=preProc)
        """
        if preProc is None:
            preProc = lambda row,fields,nameToInd: ( row, )
        if hdrPreProc is None:
            hdrPreProc = lambda row: row
        if defField is None:
            defField = copy(self.defField)
        if isinstance(csvFile,str):
            closeCsv = True
            csvFileInp = openCompressed(csvFile,"r")
            csvFile = csv.reader(csvFileInp,dialect=dialect,**dialect_options)
        else:
            closeCsv = False
            #Here we need to figure out if csvFile is just a file stream, or
            #a CSV reader already. Unfortunately, csv module does not specify
            #any common base class for CSV readers, so we have to rely on
            #tests for attribute presence
            if not (hasattr(csvFile,"dialect") and hasattr(csvFile,"line_num")):
                #this is NOT a result of calling csv.reader() or compatible interface,
                #so we assume it to be a file stream object, and call csv.reader
                #on it to create a CSV reader
                csvFile = csv.reader(csvFile,dialect=dialect,**dialect_options)

        if not hasHeader:
            firstRow = csvFile.next()
            nFields = len(firstRow)
        else:
            firstRow = None
            nFields = 0
            hdr = hdrPreProc(csvFile.next())
            #now fieldsMap is assumed to map names to SqlField, and we convert it
            #to positional map, checking first that all mapped fields are present
            #in the header
            flds_hdr = [ f.strip() for f in hdr ]
            flds_hdr_set = set(flds_hdr)
            for fldname in fieldsMap:
                if not fldname in flds_hdr_set:
                    raise ValueError,"fieldsMap argument contains field name that is "+\
                            "not found in the CSV header: %s" % (fldname,)
            fieldsMapPos = {}
            for (ifld,fldname) in enumerate(flds_hdr):
                if fldname in fieldsMap:
                    fieldsMapPos[ifld] = fieldsMap[fldname]
                    #if not already set, the field name is taken from the header
                    if fieldsMapPos[ifld].name is None:
                        fieldsMapPos[ifld].name = fldname
                else:
                    flddef = copy(defField)
                    flddef.name = fldname
                    fieldsMapPos[ifld] = flddef
            fieldsMap = fieldsMapPos
            # with header, all fields get defined in fieldsMap
            defField = None
        fields = buildSqlFields(fieldsMap=fieldsMap,nFields=nFields,defField=defField)
        table = SqlTable(name=name,fields=fields)
        self.ddl(table.createSql(),dropList=["table "+table.name])
        inserter = self.makeBulkInserterFile(table=table)
        nameToFieldInd = dict( ( (field.name,iField) for (iField,field) \
                in enumerate(fields) ) )
        if firstRow is not None:
            newRows = preProc(firstRow,fields,nameToFieldInd)
            for newRow in newRows:
                inserter(newRow)
        for row in csvFile:
            newRows = preProc(row,fields,nameToFieldInd)
            for newRow in newRows:
                inserter(newRow)
        inserter.close()
        if indices is not None and len(indices) > 0:
            self.createIndices(table=name,**indices)
            self.ddl("analyze table %s" % name,ifDialect="mysql")
        if closeCsv:
            csvFileInp.close()

    def selectAsArray(self,sql):
        """Execute SQL select and return the entire result set as Numpy record array."""
        reader = self.makeBulkReader(sql=sql)
        ret = reader.allAsArray()
        reader.close()
        return ret

    def exportAsCsv(self,sql,out,withHeader=True,bufLen=100000,
            dialect="excel-tab",
            dialect_options={"lineterminator":"\n"},
            comment=None,
            sqlAsComment=False,
            commentEscape='#',
            epilog=None):
        """Excecute SQL and export the result as CSV file.
        @param sql SQL select statement to export results of
        @param out Either file name, or file stream object, or CSV writer object
        @param withHeader If True, write the field names as the header
        @param bufLen Size (in number of records) of the internal memory buffer
        used when moving SQL result set into the output file
        @param dialect Dialect string defined in csv module
        @param dialect_options Passed to csv.writer(**dialect_options)
        @param comment if not None, this string will be printed at the top
        @param sqlAsComment if True, will print sql statement as an extra comment line
        @param commentEscape this string will be inserted at the start of every
        comment line
        @note To output any comments, out should not be a csv.writer instance
        @note We set the default lineterminator to Linux style '\n', as opposed to 
        Python's default of Windows style '\r\n'"""
        if isinstance(out,str):
            out = openCompressed(out,'w')
            doClose=True
            w = csv.writer(out,dialect=dialect,**dialect_options)
        else:
            doClose=False
            if not (hasattr(out,"dialect") and hasattr(out,"writerows")):
                #this is NOT a result of calling csv.writer() or compatible interface,
                #so we assume it to be a file stream object, and call csv.writer()
                #on it to create a CSV writer
                w = csv.writer(out,dialect=dialect,**dialect_options)
            else:
                if comment is not None or sqlAsComment:
                    raise ValueError("Illegal to write comment lines into csv.writer")
                w = out
        reader = self.makeBulkReader(sql=sql,format="list",bufLen=bufLen)
        chunks = reader.chunks()
        names = reader.fieldNames()
        if sqlAsComment:
            if comment is None:
                comment = sql
            else:
                comment = comment + '\n' + sql
        if comment is not None:
            for line in ( (commentEscape + l + '\n') for l in comment.split('\n') ):
                out.write(line)
        if withHeader:
            w.writerow(names)
        for rows in chunks:
            w.writerows(rows)
        reader.close()
        if epilog is not None:
            out.write(epilog)
        if doClose:
            out.close()

    def exportAsPivotCsv(self,sql,
            rowField,colField,valField,
            out,
            withHeader=True,
            rowFieldOut=None,
            restval=0,
            bufLen=100000,
            dialect="excel-tab",
            dialect_options={"lineterminator":"\n"},
            comment=None,
            sqlAsComment=False,
            commentEscape='#',
            epilog=None):
        """Excecute SQL and export the result as CSV file.
        @param sql SQL select statement to export results of
        @param out Either file name, or file stream object, or CSV writer object
        @param restval What to write for missing values in each cell (default is 0)
        @param withHeader If True, write the field names as the header
        @param bufLen Size (in number of records) of the internal memory buffer
        used when moving SQL result set into the output file
        @param dialect Dialect string defined in csv module
        @param dialect_options Passed to csv.writer(**dialect_options)
        @param comment if not None, this string will be printed at the top
        @param sqlAsComment if True, will print sql statement as an extra comment line
        @param commentEscape this string will be inserted at the start of every
        comment line
        @note To output any comments, out should not be a csv.writer instance
        @note We set the default lineterminator to Linux style '\n', as opposed to 
        Python's default of Windows style '\r\n'"""
        cols = sorted(self.selectAs1Col("""
        select distinct %s from 
        ( %s ) a""" % (colField,sql)))
        if rowFieldOut is None:
            rowFieldOut = rowField
        assert rowFieldOut.strip().lower() not in [ c.strip().lower() for c in cols ],\
                "%s name for row ID header name conflicts with one of the pivot column names" % (rowFieldOut,)
        names = [ rowFieldOut ] + cols
        if isinstance(out,str):
            out = openCompressed(out,'w')
            doClose=True
            w = csv.DictWriter(out, fieldnames=names, restval=restval,dialect=dialect,**dialect_options)
        else:
            doClose=False
            if not (hasattr(out,"dialect") and hasattr(out,"writerows")):
                #this is NOT a result of calling csv.writer() or compatible interface,
                #so we assume it to be a file stream object, and call csv.writer()
                #on it to create a CSV writer
                w = csv.DictWriter(out, fieldnames=names, restval=restval,dialect=dialect,**dialect_options)
            else:
                if comment is not None or sqlAsComment:
                    raise ValueError("Illegal to write comment lines into csv.writer")
                w = out

        if sqlAsComment:
            if comment is None:
                comment = sql
            else:
                comment = comment + '\n' + sql
        if comment is not None:
            for line in ( (commentEscape + l + '\n') for l in comment.split('\n') ):
                out.write(line)
        if withHeader:
            w.writerow(dict([(name,name) for name in names]))
        reader = self.makeBulkReader(sql=sql,format="list",bufLen=bufLen)
        rdrCols = reader.fieldNames(format="dict")
        # convert KeyError exception into more descriptive messages because naming mismatch between SQL
        # output fields and requested key-value fields can be a fairly typical user error
        try:
            rowFieldInd = rdrCols[rowField]
        except KeyError:
            raise ValueError("Row key field name is not found in SQL cursor recordset: %s" % (rowField,))
        try:
            colFieldInd = rdrCols[colField]
        except KeyError:
            raise ValueError("Column key field name is not found in SQL cursor recordset: %s" % (colField,))
        try:
            valFieldInd = rdrCols[valField]
        except KeyError:
            raise ValueError("Value field name is not found in SQL cursor recordset: %s" % (valField,))
        for rowKey,rows in it.groupby(reader.rows(),lambda r: r[rowFieldInd]):
            rowOut = dict(( (row[colFieldInd],row[valFieldInd]) for row in rows ))
            rowOut[rowFieldOut] = rowKey
            w.writerow(rowOut)
        reader.close()
        if epilog is not None:
            out.write(epilog)
        if doClose:
            out.close()


def sqlInList(l):
    """Create a string that can be used after SQL 'IN' keyword from a given Python sequence"""
    return "("+','.join(["%s" % x for x in l])+")"



def buildSqlFields(fieldsMap={},nFields=0,defField=None):
    """Build list of SqlField objects using field_index-to-SqlField dict or default values otherwise.
    This allows to express this pattern: I want describe a table with 10 fields, where the first field
    has name "key" and type "integer", the third field has name "value" and type varchar(10),
    and the remaining fields should be named "field_xxx" where "xxx" is a running index, and have
    type "char(50)".
    The primary use case is to create a table from CSV file where CSV file does not have the header, and
    we will only manipulate the data using a few specific fields and care to give them names and types,
    but want to carry others along as the data payload.
    @param fieldsMap Dict ( field index -> SqlField instance ) for some set of field index values
    @param nFields The total number of fields in the table (will be ajusted to the max field index 
    in fieldsMap plus one)
    @param defField SqlField instance that will be used to generate SqlField objects for field positions
    not present in fieldsMap. defField.name will be used as a prefix for new fields. If defField is None,
    SqlField(name="fld",type=char(40)) will be used.
    """
    nFields = max(list(fieldsMap.keys())+[nFields-1]) + 1
    if defField is None:
        defField = SqlField(name="fld",type="char(40)")
    fields = [] 
    for ifld in xrange(nFields): 
        fld = copy(defField)
        fld.name = "%s_%00i" % (fld.name,ifld)
        fields.append(fld)
    for (ifld,fld) in fieldsMap.items():
        fields[ifld] = fld
    return fields




class DbSqlLite(DbSql):
    """Derivative of DbSql specific for SQLite DB engine"""
    
    defField = SqlField(name="fld",type="text")
    
    def __init__(self,dbpath,strType=str,dryRun=False):
        #from pysqlite2 import dbapi2 as dbmod
        #it is a standard Python module since python 2.5:
        import sqlite3 as dbmod
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
            delattr(self,'con')

    def makeBulkInserterFile(self,**kw):
        #@todo There is a lot of room for optimizing bulk insertion
        #in SQLite. E.g. this
        #$echo "create table mytable ( col1 int, col2 int);" | sqlite3 foo.sqlite
        #$echo ".import demotab.txt mytable"  | sqlite3 foo.sqlite
        #or other methods that use true sql insert but wrap many of them in a 
        #single transction as well as play with memory settings, described e.g.
        #here: http://stackoverflow.com/questions/364017/faster-bulk-inserts-in-sqlite3
        k = {}
        if kw.has_key("bufLen"):
            k["bufLen"] = kw["bufLen"]
        tbl = kw["table"]
        if not hasattr(tbl,"insertSql"):
            tbl = SqlTable.fromDb(db=self,name=tbl)
        k["sql"] = tbl.insertSql()
        return self.makeBulkInserter(**k)
   


class DbSqlMy(DbSql):
    """Derivative of DbSql specific for MySQL DB engine"""
    
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
        self.con = self.dbmod.connect(port=kw.get("port",13306),
                                      #unix_socket=kw.get("unix_socket","/tmp/atovtchi.mysql.sock"),
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


    def makeBulkInserterFile(self,**kw):
        kw = copy(kw)
        if hasattr(kw["table"],"name"):
            kw["table"] = kw["table"].name
        return BulkInserterFileMy(self,**kw)

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
    """Derivative of DbSql specific for MonetDB DB engine"""
    
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
        """This columnar RDBMS has not explicit indices"""
        pass

    def createTableAs(self,name,select,indices=None):
        """Specialization.
        MonetDB requires 'with data' suffix"""

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
        if indices is not None and len(indices) > 0:
            self.createIndices(table=name,**indices)
            self.ddl("analyze table %s" % name,ifDialect="mysql")

def createDbSql(**kw):
    #db = DbSqlLite("/export/atovtchi/test_seq.db")
    #db = DbSqlMonet()
    db = DbSqlMy(**kw)
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
    
    def __init__(self,db,sql,bufLen=100000):
        self.db = db
        self.bufLen = bufLen
        self.sql = sql
        self.buf = []
        
    def __call__(self,record):
        self.buf.append(record)
        if len(self.buf) >= self.bufLen:
            self.flush()

    def __del__(self):
        self.db = None
            
    def flush(self):
        if len(self.buf) > 0:
            self.db.executemany(self.sql,self.buf)
            self.buf = []

    def close(self):
        self.flush()

class BulkInserterFile(object):
    
    def __init__(self,db,table,bufLen=500000,workDir="."):
        """Bulk loading from files is not part of SQL standard.
        However, every DBMS with a server process seems to have a way to do it, and it is very fast
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


def _debugTaxaCount(db):
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
            if fld_t is None:
                npy_t = 'O'
            elif fld_t == dbapi.NUMBER:
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

    def rows(self):
        """Iterate row-by-row"""
        for chunk in self.chunks():
            for row in chunk:
                yield row

    def allAsArray(self):
        """Return all (remaining) records as one numpy record array or nested list."""
        rows = self.curs.fetchall()
        if self.format == "array":
            if len(rows):
                return numpy.rec.fromrecords(list(rows),dtype=self.dt)
            else:
                return numpy.empty(0,dtype=self.dt)
        else:
            return rows

    def fieldNames(self,format="list"):
        """Return field names.
        @param format If 'list' [ default ] - return a list of names,
        if 'dict' - return a dict that maps field names to positions
        """
        ret = [ fld[0] for fld in self.curs.description ]
        if format == "list":
            return ret
        elif format == "dict":
            return dict([ (x[1],x[0]) for x in enumerate(ret) ])
        else:
            raise ValueError("Unknown value for 'format': %s" % (format,))

    def close(self):
        self.curs.close()
        self.db = None


