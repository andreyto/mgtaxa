"""Random access database for collected taxonomically characterized sequence data.
"""

from MGT.Common import *
from MGT.Taxa import *
from MGT.Sql import *
from MGT.BlastDb import BlastDb

import tables as pt

class HdfSeqInd(pt.IsDescription):
    id        = pt.Int64Col(pos=1) # primary key id, currently NCBI GI
    # [begin,begin + size) form right-open-ended range, same as python range() or STL begin(),end()
    begin     = pt.Int64Col(pos=2) # index of first element within 'seq' dataset
    size       = pt.Int64Col(pos=3) # size of sequences

class HdfSeqLoader(Options):

    def __init__(self,hdfFile,hdfGroupName,seqSizeEstimate):
        Options.__init__(self)
        self.hdfFile = hdfFile
        hdfGroup = hdfFile.createGroup(hdfFile.root,hdfGroupName)
        atomSeq = pt.StringAtom(itemsize=1)
        # Use ``a`` as the object type for the enlargeable array.
        seq = hdfFile.createEArray(hdfGroup,
            'seq',
            atomSeq,
            (0,),
            "Sequence",
            expectedrows=seqSizeEstimate,
            filters=pt.Filters(complevel=2, complib='zlib',shuffle=False)
            )
        #size of internal memory buffer of this leaf
        seq.nrowsinbuf = 1024**2 #1M
        self.seq = seq
        ind = hdfFile.createTable(hdfGroup, 'ind', HdfSeqInd, "Sequence Index")
        #size of internal memory buffer of this leaf
        ind.nrowsinbuf = 1024**2 #1M
        self.ind = ind

    def loadBlastDb(self,giFile):
        blastAlias = self.blastSelAlias
        blastDb = BlastDb()
        blastDb.makeDbAlias(blastAlias,giFile)
        #defLineTargetOnly is True because we should have already removed duplicate
        #records while collecting the sequence headers
        fastaInp = blastDb.fastaReader(dbName=blastAlias,giFile=giFile,defLineTargetOnly=True)
        indBegin = 0
        chunkSize = 1024**2
        iRec = 0
        for rec in fastaInp.records():
            gi = int(rec.getNCBI_Id())
            nSeq = 0
            indRec = self.ind.row
            for chunk in rec.seqArrays(chunkSize=chunkSize):
                nSeq += len(chunk)
                self.seq.append(chunk)
            indRec['id'] = gi
            indRec['begin'] = indBegin
            indRec['size'] = nSeq
            indRec.append()
            indBegin += nSeq
            iRec += 1
            if iRec % 10000 == 0:
                print "Done %s records of total length %s" % (iRec,indBegin)
        fastaInp.close()

    def close(self):
        self.seq.close()
        self.seq = None
        self.ind.close()
        self.ind = None
        self.hdfFile.flush()
        self.hdfFile = None

class HdfSeqReader(Options):

    def __init__(self,hdfFile,hdfGroupName):
        Options.__init__(self)
        self.hdfFile = hdfFile
        seq = hdfFile.getNode("/"+hdfGroupName+"/seq")
        seq.nrowsinbuf = 1024**2 #1M
        self.seq = seq
        ind = hdfFile.getNode("/"+hdfGroupName+"/ind")
        ind.nrowsinbuf = 1024**2 #1M
        self.ind = ind

class HdfSeqReaderSql(HdfSeqReader):

    def __init__(self,db,*l,**kw):
        HdfSeqReader.__init__(self,*l,**kw)
        self.db = db

    def loadIndToSql(self):
        db = self.db
        db.ddl("create table seq_ind (id bigint, begin bigint, size bigint)",
            dropList=["table seq_ind"])
        inserter = db.makeBulkInserterFile(table="seq_ind",bufLen=500000,workDir=self.tmpDir)
        for row in self.ind:
            inserter(row[:])
        inserter.flush()
        db.ddl("ANALYZE TABLE seq_ind",ifDialect="mysql")
        db.createIndices(names=["begin"],table="seq_ind",primary="id")
        db.ddl("ANALYZE TABLE seq_ind",ifDialect="mysql")
