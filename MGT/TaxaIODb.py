### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.TaxaTree import TaxaNode
from MGT.Sql import *



def sqlIsSubTree(aliasSub,aliasSup,withEquality=False):
    """SQL 'macro'. Return part of the WHERE clause in round braces: require that node from 'aliasSub' is in sub-tree of the node from 'aliasSup'.
    It uses 'nested sets' index fields 'lnest' and 'rnest'.
    @param aliasSub SQL select alias for the table with sub-nodes
    @param aliasSup SQL select alias for the table with super-nodes
    @param withEquality if True, a node is considered to be in the subtree of itself"""
    if withEquality:
        opmore = ">="
        opless = "<="
    else:
        opmore = ">"
        opless = "<"
    return "( %(sub)s.lnest %(opmore)s %(sup)s.lnest and %(sub)s.rnest %(opless)s %(sup)s.rnest )" % \
        dict(sub=aliasSub,sup=aliasSup,opmore=opmore,opless=opless)

class NodeStorageDb:
    """Taxonomy nodes storage in SQL DB"""
    ##@todo Implement saving and loading of 'merged' node map.
    ## Currently, we consider taxa names immutable and originating
    ## from the initial NCBI taxonomy dump file, so load them from a fixed
    ## DB table named by 'tblNames'. This assumes that all other trees that
    ## we construct are reductions of that original NCBI tree. This might change
    ## in the future, in which case we will use specific name table for each tree.
    ## It will still make sense to keep the names in a separate table because they
    ## are long and rarely needed strings.
    
    tblNames = "taxa_names"

    tblPrefix = options.taxaTreeTablePrefix

    def __init__(self,db,tableSfx):

        self.db = db
        self.tableSfx = tableSfx
        self.tblNodes = self.tblPrefix+tableSfx

    def load(self):
        db = self.db
        reader = db.makeBulkReader(sql="select * from %s" % (self.tblNodes,),bufLen=100000)
        nodes = {}
        for chunk in reader.chunks():
            for rec in chunk:
                node = TaxaNode()
                node.id = rec['id']
                node.idpar = rec['idpar']
                node.rank = rec['rank']
                node.lnest = rec['lnest']
                node.rnest = rec['rnest']
                node.depth = rec['depth']
                node.seq_len = rec['seq_len']
                node.seq_len_tot = rec['seq_len_tot']
                node.idlevel = rec['idlevel']
                node.divid = rec['divid']
                nodes[node.id] = node
        reader.close()
        return dict(nodes=nodes,merged={})

    def save(self,tree):
        db = self.db
        db.ddl("""
        create table %s
        (
        id integer,
        idpar integer,
        lnest integer,
        rnest integer,
        seq_len bigint,
        seq_len_tot bigint,
        idlevel tinyint,
        depth tinyint,
        rank char(20),
        divid tinyint
        )
        """ % (self.tblNodes,),
        dropList=["table %s" % (self.tblNodes,)])

        inserter = db.makeBulkInserterFile(table=self.tblNodes,bufLen=500000)
        for n in tree.iterDepthTop():
            inserter((n.id,
                      n.idpar,
                      n.lnest,
                      n.rnest,
                      n.seq_len if hasattr(self,"seq_len") else None,
                      n.seq_len_tot if hasattr(self,"seq_len") else None,
                      n.idlevel,
                      n.depth,
                      n.rank,
                      n.divid))
        inserter.flush()
        db.ddl("ANALYZE TABLE %s" % (self.tblNodes,),ifDialect="mysql")
        db.createIndices(table=self.tblNodes,primary="id",
                names=["idpar","lnest","rnest","idlevel","depth","rank","divid"],
                attrib={"lnest":{"unique":True},"rnest":{"unique":True}})

    def saveName(self,tree):
        """Save node names into a separate table called by self.tblNames attribute.
        @param Tree object"""
        #this method could be also implemented by a single call to self.saveAttributes()
        db = self.db
        db.ddl("""\
        create table %s
        (
        id integer,
        name  varchar(180)
        )
        """ % (self.tblNames,),
        dropList=["table %s" % (self.tblNames,)])
        inserter = db.makeBulkInserterFile(table=self.tblNames,bufLen=500000)
        for n in tree.iterDepthTop():
            inserter((n.id,
                      n.name))
        inserter.flush()
        db.ddl("ANALYZE TABLE %s" % (self.tblNames,),ifDialect="mysql")
        db.createIndices(table=self.tblNames,primary="id")

    def loadAttribute(self,tree,name,sql,default=None,setDefault=True,ignoreKeyError=True,typeCast=None):
        """Execute 'sql' which should return a unique mapping taxid -> value and assign result to each node.
        @param name - name of the new tree node attribbute to set
        @param sql - statement to execute. It must return (id,value) pairs (in that order) with 'id' corresponding to tree node id's.
        Actual column names do not matter.
        @param default - assign this value to those nodes for which taxid is not present in the 'sql' result set (if 'setDeafult' is True)
        @param setDefault - if False, do not assign default value to nodes
        @param ignoreKeyError - if True, do not raise exception if 'sql' results set contains an id that is not present in the tree.
        """
        db = self.db
        if setDefault:
            tree.setAttribute(name,default)
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
                    setattr(tree.getNode(int(row[0])),name,typeCast(row[1]))
                except KeyError:
                    if not ignoreKeyError:
                        raise
        
    
    def saveAttributes(self,tree,nameToSql,table):
        """Save attribute values of tree nodes as a new table indexed by node id.
        @param tree - tree to work on
        @param nameToSql - mapping of attribute names into SQL, by example: {'rankUpper':{'name':'rank_upper','type':'char(20)','createIndex':True,default=0}}. If 'name' is not present, the key of the dictionary will be used as is for SQL column
        name. That might lead to problems if SQL server converts all names to lower case, for instance. It is not a good
        idea to rely on identifier case in SQL code anyway.
        @param table - SQL table name to create
        """
        db = self.db
        #protect against typos in keyword names
        attribKeywords = set(('name','type','createIndex','default'))
        for name in nameToSql:
            for k in nameToSql[name].keys():
                if k not in attribKeywords:
                    raise ValueError("Unknown attribute keyword: " + name + "[" + k + "]")
        sqlColDef = ",\n".join([ "%(name)s %(type)s" % {'name':v.get('name',k),'type':v['type']} for (k,v) in nameToSql.items() ])
        db.ddl("""
        create table %s
        (
        id integer,
        %s
        )
        """ % (table,sqlColDef),
        dropList=["table %s" % (table,)])

        inserter = db.makeBulkInserterFile(table=table,bufLen=500000)
        #python guarantees the same order as from .values()
        namesAndDefaults = [ (name,val.get('default',None)) for (name,val) in nameToSql.items() ]
        for n in tree.iterDepthTop():
            vals = [ getattr(n,name,default) for (name,default) in namesAndDefaults ]
            # If all values are NULL (None), do not insert the record.
            # Optimized for common case of all values non-NULL.
            if None in vals:
                allNone = True
                for v in vals:
                    if v is not None:
                        allNone = False
                        break
                if allNone:
                    continue
            vals.insert(0,n.id)
            inserter(vals)
        inserter.flush()
        namesToIndex = [ nameToSql[name].get('name',name) for name in nameToSql if nameToSql[name].get('createIndex',False) ]
        db.createIndices(table=table,primary="id",names=namesToIndex)


    def loadSeqLen(self,
            tree,
            name="seq_len",
            sql="select taxid,seq_len from taxa_seq_len"):
        """Load sequence length from DB and also assign accumulated subtree sequence to attribute name+'_tot'
        Parameters have the same meaning as for for loadAttribute().
        Default value of 0L is used."""
        self.loadAttribute(tree=tree,
                           name=name,
                           sql=sql,
                           default=0L,
                           setDefault=True,
                           ignoreKeyError=True,
                           typeCast=long)
        tree.setTotal(srcAttr=name,dstAttr=name+"_tot")

class LevelsStorageDb:
    """Taxonomy levels storage in SQL DB"""

    tblPrefix = options.taxaLevelsTablePrefix

    def __init__(self,db,tableSfx=""):

        self.db = db
        self.tableSfx = tableSfx
        self.tblLevels = self.tblPrefix+tableSfx

    def save(self,taxaLevels):
        db = self.db
        db.ddl("""
        create table %s
        (
        id tinyint,
        level char(20),
        is_linn bool,
        pos tinyint
        )
        """ % (self.tblLevels,),
        dropList=["table %s" % (self.tblLevels,)])
        levelIds = taxaLevels.getIdToLevel() # has no_rank
        levelPos = taxaLevels.getLevelsPos(order="ascend") # does not have no_rank
        inserter = db.makeBulkInserter(table=self.tblLevels,bufLen=500000)
        for (idlevel,level) in sorted(levelIds.items()):
            is_linn = level in levelPos
            inserter((idlevel,
                      level,
                      is_linn,
                      levelPos[level] if is_linn else None))
        inserter.flush()
        db.ddl("ANALYZE TABLE %s" % (self.tblLevels,),ifDialect="mysql")
        db.createIndices(table=self.tblLevels,primary="id",
                names=["level","is_linn","pos"],
                attrib={"level":{"unique":True},"pos":{"unique":True}})

