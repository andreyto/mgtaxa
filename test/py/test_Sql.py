### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Sql import *

db = DbSqlMy()

db.ddl("""
    create table tmp_test_types
        (f_bool bool, f_tinyint tinyint, f_med mediumint, f_int int,
        f_big bigint, f_char char(5), f_float float)
    """,
    dropList=["table tmp_test_types"])

db.ddl("""
    insert into tmp_test_types
        values
        (%s)
    """ % (','.join(['1']*7),))
db.ddl("""
    insert into tmp_test_types
        values
        (%s)
    """ % (','.join(['2']*7),))

expect = numpy.rec.fromrecords([(1, 1, 1, 1, 1, '1', 1.0), (2, 2, 2, 2, 2, '2', 2.0)],
    dtype=[('f_bool', '|i1'), ('f_tinyint', '|i1'),
           ('f_med', '<i4'), ('f_int', '<i4'),
           ('f_big', '<i8'), ('f_char', '|S5'),
           ('f_float', '<f8')])

reader = db.makeBulkReader(sql="select * from tmp_test_types",bufLen=100)

for chunk in reader.chunks():
    print chunk
    print chunk.shape, chunk.dtype
    print reader.nrowsFetched()
    assert numpy.all(chunk == expect)
