"""Classes that specialize Taxonomy trees to be represented in SQL database and operations on them."""

from MGT.Common import *
from MGT.TaxaTree import *
from MGT.Sql import *

class TaxaTreeDb(TaxaTree):
    """Taxonomy tree that generates and loads data for/from SQL DB."""
    
    def __init__(self,storage,db,tableSfx=None):
        """Initialize new TaxaTreeDb instance.
        @param db - DbSql object. 
        @param tablePrefix - DB tables with data for this tree will have this prefix"""
        TaxaTree.__init__(self,storage=storage)
        self.db = db
        self.tableSfx = tableSfx

    def loadAttribute(self,name,sql,default=None,setDefault=True,ignoreKeyError=True,typeCast=None):
        """Execute 'sql' which should return a unique mapping taxid -> value and assign result to each node.
        @param name - name of the new tree node attribute to set
        @param sql - statement to execute. It must return (id,value) pairs (in that order) with 'id' corresponding to tree node id's.
        Actual column names do not matter.
        @param default - assign this value to those nodes for which taxid is not present in the 'sql' result set (if 'setDeafult' is True)
        @param setDefault - if False, do not assign default value to nodes
        @param ignoreKeyError - if True, do not raise exception if 'sql' results set contains an id that is not present in the tree.
        """
        db = self.db
        if setDefault:
            self.setAttribute(name,default)
        if typeCast is None:
            typeCast = lambda v: v
        curs = db.execute(sql)
        curs.arraysize = 100000
        while True:
            rows = curs.fetchmany()
            if not rows:
                break
            for row in rows:
                try:
                    setattr(self.getNode(int(row[0])),name,typeCast(row[1]))
                except KeyError:
                    if not ignoreKeyError:
                        raise
        
    
    def loadSeqLen(self,name="seq_len",
                         sql="select taxid,seq_len from taxa_seq_len")
        """Load sequence length from DB and also assign accumulated subtree sequence to attribute name+'_tot'
        Parameters have the same meaning as for for loadAttribute().
        Default value of 0L is used."""
        self.loadAttribute(name=name,
                           sql=sql,
                           default=0L,
                           setDefault=True,
                           ignoreKeyError=True,
                           typeCast=long):
        self.setTotal(srcAttr=name,dstAttr=name+"_tot")


class TaxaLevelsDb(TaxaLevels,Options):
    """Load our selection of taxonomic ranks into SQL database, create column and row aggregates.
    """

    def __init__(self,db,taxaTree):
        Options.__init__(self)
        self.taxaTree = taxaTree
        TaxaLevels.__init__(self)
        if self.taxaTree is not None:
            self.setLevels(taxaTree)
        self.db = db

    def getLevelColumns(self,order="ascend"):
        return [ "ti_"+name for name in self.getLevelNames(order=order) ]

    def getLevelColumnsComma(self,alias=None,order="ascend"):
        """Return a string with a comma separated list of level column names.
        @param alias - If not None, should be a one-letter SQL table alias to prepend each column name with, e.g. 'b.ti_family'
        @param order - return in the order of ascending or descending ranks."""
        cols = self.getLevelColumns(order=order)
        if alias is not None:
            alias = alias.strip()
            assert len(alias) == 1 and alias.isalpha()
            cols = [ alias+'.'+c for c in cols ]
        return ','.join(cols)

    def loadTaxLevelsRows(self):
        """@todo create table taxa_tree = taxa_node + unclass + dist_top + seq_len(node+subtree)
        seq_len calc on memory tree after loading leaf seq_len from all_src"""
        self.db.ddl("""\
        create table taxa_level_row
        (
        taxid integer,
        unclass bool,
        idlevel tinyint,
        partaxid integer,
        distance tinyint
        )
        """,
        dropList=["table taxa_level_row"])
        
        inserter = self.db.makeBulkInserterFile(table="taxa_level_row",bufLen=500000,workDir=self.tmpDir)
        taxaTree = self.taxaTree
        for node in taxaTree.iterDepthTop():
            for (idlevel,partaxid,distance) in self.lineageKeys(node,getDistance=True):
                inserter((node.id,int(node.isUnclassified()),idlevel,partaxid,distance))
        inserter.flush()
        self.db.createIndices(table="taxa_level_row",
            names=["taxid","unclass","partaxid","idlevel","distance"],
            primary="taxid,idlevel,partaxid")

    def loadTaxLevelsColumns(self):
        taxaTree = self.taxaTree
        ti_cols = self.getLevelColumns()
        self.db.ddl("""\
        create table taxa_level_col
        (
        taxid integer,
        %s
        )
        """ % (",\n".join(["%s integer" % (col,) for col in ti_cols]),),
        dropList=["table taxa_level_col"])

        inserter = self.db.makeBulkInserterFile(table="taxa_level_col",bufLen=500000,workDir=self.tmpDir)
        for node in taxaTree.iterDepthTop():
            inserter([node.id]+self.lineageFixedList(node))
        inserter.flush()
        self.db.createIndices(table="taxa_level_col",
            names=ti_cols,
            primary="taxid")

    def makeStatsTables(self):
        db = self.db
        
        db.createTableAs("taxa_level_row_len","""
        select a.*,b.seq_len from taxa_level_row a,taxa_seq_len b where a.taxid = b.taxid
        """)

        db.createIndices(names=["taxid","unclass","partaxid","idlevel","distance"],table="taxa_level_row_len")

        db.executeAndPrint("""
            select
            idlevel,count(*) as cnt,
            avg(seq_len) as seq_len
            from (
                select idlevel,partaxid,count(*) as cnt,sum(seq_len) as seq_len
                from taxa_level_row_len
                group by idlevel,partaxid) a
            group by idlevel
        """)

        ti_cols = self.getLevelColumns()
        ti_group = list(ti_cols)
        ti_group.reverse()
        ti_group_comma = ','.join(ti_group)
        group_names = self.getLevelNames()
        group_names.reverse()
        
        db.createTableAs("taxa_level_col_len","""
        select a.*,b.seq_len from taxa_level_col a,taxa_seq_len b where a.taxid = b.taxid
        """)

        db.createIndices(table="taxa_level_col_len",names=ti_cols,primary="taxid")

        ## We replace NULL of taxid values with 0, so that
        ## the MySQL OLAP modifier 'with rollup' would work correctly
        ## (it uses NULL to mark rows with totals)
        
        db.createTableAs("taxa_level_col_gr","""
        select %s,sum(seq_len) as seq_len,count(*) as taxid_cnt
        from taxa_level_col_len a
        group by %s
        """ % (",".join(["COALESCE(a.%s,0) AS %s" % (ti_gr,ti_gr) for ti_gr in ti_group]) ,",".join(["a."+ti_gr for ti_gr in ti_group])))
        
        db.createIndices(table="taxa_level_col_gr",names=ti_cols)


        db.createTableAs("taxa_level_col_rep_1","""
        select %s,sum(seq_len) as seq_len,sum(taxid_cnt) as taxid_cnt
        from taxa_level_col_gr
        group by %s with rollup
        """ % (ti_group_comma,ti_group_comma))
        
        db.ddl("""ALTER TABLE taxa_level_col_rep_1 ADD id INTEGER auto_increment  PRIMARY KEY""")
        #db.createIndices(table="taxa_level_col_rep_1",names=ti_cols)

        ## Add columns with string names for each ti_xxx column through 'left join'

        from string import ascii_lowercase
        alias_a = ascii_lowercase[0]
        alias_joins = ascii_lowercase[1:]
        assert len(ti_group) <= len(alias_joins)
        joins = "\n".join([ "LEFT JOIN taxa_names %s ON %s.%s = %s.taxid" % (ali_join,alias_a,ti_col,ali_join)
                        for (ali_join,ti_col) in zip(alias_joins,ti_group) ])

        cols = ",\n".join(["%s,%s.name AS nm_%s" % (ti_col,ali_join,gr_name) for (ti_col,ali_join,gr_name) in zip(ti_group,alias_joins,group_names)])
                    
        db.createTableAs("taxa_level_col_rep","""
        SELECT
        id,
        %s,
        seq_len,taxid_cnt
        FROM taxa_level_col_rep_1 %s
        %s
        ORDER BY id
        """ % (cols,alias_a,joins))

        db.createIndices(table="taxa_level_col_rep",names=ti_group+["id"])
        
        db.executeAndPrint("""
        select ti_superkingdom,nm_superkingdom,seq_len,taxid_cnt
        from taxa_level_col_rep
        where ti_phylum is NULL
        """)
    
        #db.execute("""
        #select *
        #from taxa_level_col_rep
        #order by id
        #into outfile 'taxa_level_col_rep.csv'
        #fields
        #terminated by '|'
        #optionally enclosed by '"'
        #""")

        db.executeAndPrint("""
        select ti_phylum,nm_phylum,ti_class,nm_class,ti_order,nm_order,seq_len,taxid_cnt
        from taxa_level_col_rep
        where ti_family is NULL and ti_order is not NULL and ti_class is not NULL
        """)
        
        db.executeAndPrint("""
        select ti_phylum,nm_phylum,ti_class,nm_class,ti_order,nm_order,ti_family,nm_family,seq_len,taxid_cnt
        from taxa_level_col_rep where ti_superkingdom = 2 and ti_genus is NULL and ti_family is not NULL and
        ti_order is not NULL and ti_class is not NULL
        """)

        db.executeAndPrint("""
        select ti_phylum,nm_phylum,ti_class,nm_class,ti_order,nm_order,ti_family,nm_family,seq_len,taxid_cnt
        from taxa_level_col_rep where ti_superkingdom = 2 and ti_genus is NULL and ti_family is not NULL and
        ti_family <> 0 and ti_order is not NULL and ti_class is not NULL
        """)

