### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import MonetSQLdb
import sys

dbh = MonetSQLdb.Connection(dbname='mgtaxa')
#if len(sys.argv) > 3:
    #dbh = MonetSQLdb.Connection(dbfarm=sys.argv[1],dbname=sys.argv[2],version=int(sys.argv[3]))
#else:
    #dbh = MonetSQLdb.Connection(port=sys.argv[1],dbname=sys.argv[2])

cursor = dbh.cursor();

try:
    cursor.execute('create table python_table (i smallint,s string);');
    #this works
    s = ((0, 'row1'), (1, 'row2'))
    x = cursor.executemany("insert into python_table VALUES (%s, %s);", s)
    print(x);
    #this does not work
    s = [[0, 'row1'], [1, 'row2']]
    x = cursor.executemany("insert into python_table VALUES (%s, %s);", s)
    print(x);
finally:
    cursor.execute('drop table python_table;');
