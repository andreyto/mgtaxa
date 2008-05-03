from MGT.Common import *
from MGT.TaxaTree import TaxaNode
from MGT.Sql import *


class NodeStorageDb:
    """Taxonomy nodes storage in SQL DB"""

    ## Currently, we consider taxa names immutable and originating
    ## from the initial NCBI taxonomy dump file, so load them from a fixed
    ## DB table named by 'tblNames'. This assumes that all other trees that
    ## we construct are reductions of that original NCBI tree. This might change
    ## in the, in which case we will use specific name table for each tree.
    ## It will still make sense to keep the names in a separate table because they
    ## are long and rarely needed strings.
    ## Currently we just expect the 'tblNames' to be loaded into the database already.
    
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
                nodes[node.id] = node
        reader.close()
        return nodes

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
        rank char(20)
        )
        """ % (self.tblNodes,),
        dropList=["table %s" % (self.tblNodes,)])

        inserter = db.makeBulkInserterFile(table=self.tblNodes,bufLen=500000)
        for n in tree.iterDepthTop():
            inserter((n.id,
                      n.idpar,
                      n.lnest,
                      n.rnest,
                      n.seq_len,
                      n.seq_len_tot,
                      n.idlevel,
                      n.depth,
                      n.rank))
        inserter.flush()
        db.ddl("ANALYZE TABLE %s" % (self.tblNodes,),ifDialect="mysql")
        db.createIndices(table=self.tblNodes,primary="id",names=["idpar","lnest","rnest","idlevel","depth","rank"])
