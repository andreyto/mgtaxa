from MGT.Sql import *

#input with headr
csvHdrInp = pjoin(options.testDataDir,"sql","apis.hdr.csv")
#input w/o header
csvNhdrInp = pjoin(options.testDataDir,"sql","apis.nhdr.csv")

#db = DbSqlLite(dbpath=":memory:")
db = DbSqlLite(dbpath="tmp.db.sqlite")

db.createTableFromCsv(name="test_csv_hdr",
        csvFile=csvHdrInp,
        fieldsMap={"f1":SqlField(name="field1",type="integer")},
        defField=None,
        hasHeader=True,
        dialect="excel-tab",
        indices=None)

db.exportAsCsv(sql="select * from test_csv_hdr",
        out="test_csv_hdr.csv",
        withHeader=True,
        bufLen=100000,
        dialect="excel-tab")

#Copy "expected" by directly opening the file with redirected output in your editor (e.g. VIM),
#because Linux terminal clipboard converts tabs to spaces.
expected = """- seqid	field1	f2	f3	f4	f5	f6	f7	f8	f9	f10	f11	f12	f13	f14	f15	f16	f17	f18	f19	f20	f21	f22	f23	f24
?        ----
+ seqid	f1	f2	f3	f4	f5	f6	f7	f8	f9	f10	f11	f12	f13	f14	f15	f16	f17	f18	f19	f20	f21	f22	f23	f24
"""
received = diffFiles("test_csv_hdr.csv",csvHdrInp,asText=True)
print "received:"
print received,
print "expected:"
print expected,
assert received == expected

db.createTableFromCsv(name="test_csv_nhdr",
        csvFile=csvNhdrInp,
        fieldsMap={1:SqlField(name="field1",type="integer")},
        defField=None,
        hasHeader=False,
        dialect="excel-tab",
        indices=None)

db.exportAsCsv(sql="select * from test_csv_nhdr",
        out="test_csv_nhdr.csv",
        withHeader=True,
        bufLen=100000,
        dialect="excel-tab")

expected = """- fld_0	field1	fld_2	fld_3	fld_4	fld_5	fld_6	fld_7	fld_8	fld_9	fld_10	fld_11	fld_12	fld_13	fld_14	fld_15	fld_16	fld_17	fld_18	fld_19	fld_20	fld_21	fld_22	fld_23	fld_24
"""
received = diffFiles("test_csv_nhdr.csv",csvNhdrInp,asText=True)
print "received:"
print received,
print "expected:"
print expected,
assert received == expected

sql = "select f7,f10,count(*) as cnt from test_csv_hdr group by f7,f10 order by f7,f10"
rowField = "f7"
colField = "f10"
valField = "cnt"

db.exportAsPivotCsv(sql=sql,
        out="test_pivot.csv",
        rowField=rowField,colField=colField,valField=valField,
        withHeader=True,
        bufLen=100000,
        dialect="excel-tab")
