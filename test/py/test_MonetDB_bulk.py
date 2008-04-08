import MonetSQLdb
import sys

bulkDat = "/home/atovtchi/work/mgtaxa/test_data/bulk.dat"

dbh = MonetSQLdb.Connection(dbname='mgtaxa')

def execute(sql):
    cursor = dbh.cursor()
    cursor.execute(sql)
    return cursor

def setup():

    try:
        execute("drop table test_bulk").close()
    except StandardError, msg:
        print "Warning: ", msg
    execute("""
        create table test_bulk
        (
        iid integer,
        gi bigint,
        taxid integer,
        src_db varchar(1),
        project varchar(4),
        seq_len bigint,
        acc varchar(20),
        kind varchar(2),
        seq_hdr varchar(40)
        )
    """).close()


def insertPy():
    for iBatch in range(10):
        #COPY 5 OFFSET 5 RECORDS INTO my_test FROM stdin USING DELIMITERS '|','\n' ;
        execute("copy %s offset %s records into test_bulk from '%s'" % (50000,iBatch*50000,bulkDat,)).close()
        cursor = execute("select count(*) from (select taxid from test_bulk group by taxid) a")
        print cursor.fetchall()
        cursor.close()


def insertPyOneCurs():
    cursor = dbh.cursor()
    for iBatch in range(10):
        #COPY 5 OFFSET 5 RECORDS INTO my_test FROM stdin USING DELIMITERS '|','\n' ;
        cursor.execute("copy %s offset %s records into test_bulk from '%s'" % (50000,iBatch*50000,bulkDat,))
        cursor.execute("select count(*) from (select taxid from test_bulk group by taxid) a")
        print cursor.fetchall()
    cursor.close()



#if len(sys.argv) > 3:
    #dbh = MonetSQLdb.Connection(dbfarm=sys.argv[1],dbname=sys.argv[2],version=int(sys.argv[3]))
#else:
    #dbh = MonetSQLdb.Connection(port=sys.argv[1],dbname=sys.argv[2])

#cursor = dbh.cursor();

setup()
insertPyOneCurs()

dbh.close()