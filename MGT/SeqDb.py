"""Random access database for sequence data.
We store all sequence data in one flat character array compressed with zlib.
A separate array is built to index the sequence data.
Each index record is a tuple (id,begin,size) and thus describes an arbitrary contiguous chunk
of sequence and its identifier.
Initially we load sequences from NCBI and the original index represents individual NCBI sequence
records.
Later we can build other indices in separate HDF files to access sequence data from 'seq' dataset
in this file. E.g. an index that describes a concatenation of concequtive original sequences
with some sequences (e.g. 16S rRNA) omitted.
"""

from MGT.Common import *
from MGT.Taxa import *
from MGT.Sql import *
from MGT.BlastDb import BlastDb

from MGT.Hdf import *

from itertools import izip

class HdfSeqInd(pt.IsDescription):
    id        = pt.Int32Col(pos=1) # primary key id, our internal
    # [begin,begin + size) form right-open-ended range, same as python range() or STL begin(),end()
    begin     = pt.Int64Col(pos=2) # index of first element within 'seq' dataset
    size       = pt.Int64Col(pos=3) # size of sequences

class HdfSeqLoader(MGTOptions):

    def __init__(self,hdfFile,hdfGroup,seqSizeEstimate,seqNumEstimate):
        MGTOptions.__init__(self)
        self.hdfFile = hdfFile
        atomSeq = pt.StringAtom(itemsize=1)
        # Use ``a`` as the object type for the enlargeable array.
        seq = hdfFile.createEArray(hdfGroup,
            'seq',
            atomSeq,
            (0,),
            "Sequence",
            expectedrows=seqSizeEstimate,
            filters=pt.Filters(complevel=2, complib='zlib',shuffle=False),
            createparents=True
            )
        #size of internal memory buffer of this leaf
        seq.nrowsinbuf = 1024**2 #1M
        self.seq = seq
        ind = hdfFile.createTable(hdfGroup, 
            'ind', 
            HdfSeqInd, 
            "Sequence Index",
            expectedrows=seqNumEstimate,
            createparents=True
            )
        #size of internal memory buffer of this leaf
        ind.nrowsinbuf = 1024**2 #1M
        self.ind = ind

    def loadBlastDb(self,giIdFile):
        """Load specified subset of BLAST DB into HDF dataset.
        @param giIdFile - file with pairs 'gi id' on each line:
        gi will be used to pull sequence from BLAST DB, id will be used as 'id' field in HDF sequence index."""
        (outGi,giFile) = makeTmpFile(mode="w",bufsize=2**20,dir=self.tmpDir,prefix="hdfSeqLoader_",suffix=".gi",createParents=True)
        giFile = os.path.abspath(giFile)
        inpGiId = open(giIdFile,"r",2**20)
        for line in inpGiId:
            outGi.write(line.split()[0]+"\n")
        inpGiId.close()
        outGi.close()
        blastAlias = self.blastSelAlias
        blastDb = BlastDb()
        blastDb.makeDbAlias(blastAlias,giFile)
        #defLineTargetOnly is True because we should have already removed duplicate
        #records while collecting the sequence headers
        fastaInp = blastDb.fastaReader(dbName=blastAlias,giFile=giFile,defLineTargetOnly=True)
        inpGiId = open(giIdFile,"r",2**20)
        indBegin = 0
        chunkSize = 1024**2
        iRec = 0
        for (rec,giIdLine) in izip(fastaInp.records(),inpGiId):
            gi = int(rec.getNCBI_Id())
            giCheck,id = giIdLine.split()
            giCheck = int(giCheck)
            id = int(id)
            assert gi == giCheck, "GI mismatch between BLAST DB FASTA stream (%i) and GI input list (%i)" % (gi,giCheck)
            nSeq = 0
            indRec = self.ind.row
            for chunk in rec.seqArrays(chunkSize=chunkSize):
                nSeq += len(chunk)
                self.seq.append(chunk)
            indRec['id'] = id
            indRec['begin'] = indBegin
            indRec['size'] = nSeq
            indRec.append()
            indBegin += nSeq
            iRec += 1
            if iRec % 10000 == 0:
                print "Done %s records of total length %s" % (iRec,indBegin)
        self.ind.flush()
        self.seq.flush()
        fastaInp.close()
        inpGiId.close()
        os.remove(giFile)

    def close(self):
        self.seq = None
        self.ind = None
        self.hdfFile.flush()
        self.hdfFile = None

class HdfSeqReader(MGTOptions):

    def __init__(self,hdfFile,hdfGroup):
        MGTOptions.__init__(self)
        self.hdfFile = hdfFile
        seq = hdfFile.getNode(hdfGroup,"seq")
        seq.nrowsinbuf = 1024**2 #1M
        self.seq = seq
        ind = hdfFile.getNode(hdfGroup,"ind")
        ind.nrowsinbuf = 1024**2 #1M
        self.ind = ind

    def close(self):
        self.seq = None
        self.ind = None
        self.hdfFile.flush()
        self.hdfFile = None
                                                        
class HdfSeqReaderSql(HdfSeqReader):

    def __init__(self,db,*l,**kw):
        HdfSeqReader.__init__(self,*l,**kw)
        self.db = db

    def loadIndToSql(self):
        db = self.db
        db.ddl("create table seq_ind (id int, begin bigint, size bigint)",
            dropList=["table seq_ind"])
        inserter = db.makeBulkInserterFile(table="seq_ind",bufLen=500000,workDir=self.tmpDir)
        for row in self.ind:
            inserter(row[:])
        inserter.flush()
        db.ddl("ANALYZE TABLE seq_ind",ifDialect="mysql")
        db.createIndices(names=["begin"],table="seq_ind",primary="id")
        db.ddl("ANALYZE TABLE seq_ind",ifDialect="mysql")


def hdfMakeActiveSeqInd(db,hdfFile,hdfPath):
    """Create new HdfSeqInd dataset in the current file that describes the "active sequence" view of HdfSeq sequence.
    This lists all sequences with the same taxid consequtively, setting 'id' to 'taxid' with the exception 
    that sequences not in 'act_seq' table are excluded."""
    where,name = hdfSplitPath(hdfPath)
    seq_cnt = long(db.selectScalar("select count(*) from act_seq"))
    ind = hdfFile.createTable(where, 
            name, 
            HdfSeqInd, 
            "Active Sequence Index",
            expectedrows=seq_cnt,
            createparents=True)
    #size of internal memory buffer of this leaf
    ind.nrowsinbuf = 1024**2 #1M

    # Our ordering relies on the precondition that seq_ind was created in tree order
    inp = db.exportToStream(\
    """SELECT  c.taxid,a.id,a.begin,a.size""",
    """
    FROM            seq_ind a,
                    act_seq b,
                    act_src c
    WHERE           b.id_src = c.id AND
                    b.id = a.id
    ORDER BY        a.begin
    """,
    fieldsTerm=' ',linesTerm=r'\n')
    iter = ( (int(rec[0]),int(rec[1]),long(rec[2]),long(rec[3])) for rec in (line.split() for line in inp) )
    indRec = ind.row
    nRec = 0
    for (taxid,id,begin,size) in iter:
        indRec['id'] = taxid
        indRec['begin'] = begin
        indRec['size'] = size
        indRec.append()
        if nRec % 10000 == 0:
            print "Indexed %s records" % nRec
        nRec += 1
    ind.flush()

