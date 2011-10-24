"""Post-prediction analysis"""

from MGT.Sql import *

#TMP:
taxaGrpRanks = linnMainRanks[1:] 

def mergeStats(samples,dbOut):
    """Merge stats SQLite tables of individual samples.
    @param samples Sequence of tuples (sample id,SQLite stats DB file),
    where sample id is arbitrary unique word that will be used as
    a field value to identify each sample in a merged DB, and
    DB file is a corresponding stats DB file created by the taxonomic
    prediction algorithm.
    @param dbOut Output DB file
    """
    db = DbSqlLite(dbpath=dbOut)
    for (iSample,sample) in enumerate(samples):
        db.ddl("attach database '%s' as dbin" % (sample[1],))
        #ranks = ("",) + linnMainRanks 
        ranks = taxaGrpRanks
        for rank in ranks:
            tableBase = "scaff_pred_filt_grp_%s" % (rank,)
            #Using main. below is critical otherwise
            #db.creatTableAs wacks the dbin. table because
            #w/o the database scope SQLite resolves table name
            #into the last attached DB.
            tableOut = "main." + tableBase
            tableInp = "dbin." + tableBase
            sqlArgs = dict(tableInp=tableInp,
                    id_samp=sample[0])
            sql = """select '%(id_samp)s' as id_samp,
            *
            from %(tableInp)s""" % sqlArgs
            if iSample == 0:
                db.createTableAs(tableOut,sql)
            else:
                db.ddl("""insert into %s
                %s""" % (tableOut,sql))

        db.ddl("detach database dbin")
    db.close()

def exportAsPivots(dbInp,csvOutBase,minWeight=0):
    db = DbSqlLite(dbpath=dbInp)
    ranks = taxaGrpRanks
    for rank in ranks:
        table = "scaff_pred_filt_grp_%s" % (rank,)
        sql = """select id_samp,clade,sum_weight 
        from %(table)s 
        where clade in 
        (
        select clade from %(table)s 
        group by clade
        having max(sum_weight) >= %(minWeight)s
        )
        order by id_samp""" % dict(table=table,minWeight=minWeight)
        rowField = "id_samp"
        colField = "clade"
        valField = "sum_weight"
        out = "%s_%s.csv" % (csvOutBase,rank)
        db.exportAsPivotCsv(sql=sql,
                out=out,
                rowField=rowField,
                colField=colField,
                valField=valField)
    db.close()

def loadEnvData(csvInp,dbOut,preProc,hdrPreProc):
    db = DbSqlLite(dbpath=dbOut)
    name = "env"
    db.createTableFromCsv(name,csvFile=csvInp,hasHeader=True,
            dialect="excel-tab",indices={"names":("id_samp",)},preProc=preProc,hdrPreProc=hdrPreProc)

def exportEnvData(dbInp,csvOutBase):
    db = DbSqlLite(dbpath=dbInp)
    ranks = taxaGrpRanks
    for rank in ranks:
        table_clades = "scaff_pred_filt_grp_%s" % (rank,)
        sql = """select * from env where id_samp in 
        ( select id_samp from %(table_clades)s )
        order by id_samp""" % dict(table_clades=table_clades)
        out = "%s_%s.csv" % (csvOutBase,rank)
        db.exportAsCsv(sql=sql,
                out=out)
    db.close()

