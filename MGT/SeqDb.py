"""Random access database for collected taxonomically characterized sequence data.
"""

from MGT.Common import *
from MGT.Taxa import *
from MGT.Sql import *
from MGT.BlastDb import BlastDb

import tables as pt

class HdfSeqInd(pt.IsDescription):
    id        = Int64Col(pos=1) # primary key id, currently NCBI GI
    # [begin,end) form open-ended range, same as python range() or STL begin(),end()
    begin     = Int64Col(pos=2) # index of first element within 'seq' dataset
    end       = Int64Col(pos=3) # index of last + 1 element within 'seq' dataset 

class HdfSeqLoader(Options):

    def __init__(self,db,hdfFile,hdfGroup,seqSizeEstimate):
        Options.__init__(self)
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
        #records when collecting the sequence headers
        fastaInp = blastDb.fastaReader(dbName=blastAlias,defLineTargetOnly=True)
        for rec in fastaInp.records():
            pass


    def mergeSelWithSeq(self,skipSeq=False):
        #from itertool import izip
        outFasta = gzip.open(self.selFastaFile,'w',compresslevel=4)
        inpDump = open(self.selDumpFile,'r')
        selGiFile = os.path.abspath(self.selGiFile)
        fldsDump = "gi,taxid,src_db,kind,project,cat,stage,src_type,iid,seq_len,divid,rank".split(',')
        pipe = Popen(("fastacmd -i %s -d %s" % (selGiFile,self.srcDbNameAlias)).split(), 
                     cwd=self.blastDataDir, env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        #inpFasta = readFastaRecords(pipe,readSeq=True)
        FGI_SKIP  = 0x01
        FGI_WRITE = 0x02
        FGI_MISM  = 0x04
        giSeen = {}
        iRec = 0
        skip = True
        mismatchRun = 0 # how many gi mismatches in a row
        for line in pipe:
            try:
                if line.startswith(">"):
                    skip = False
                    #header line can be:
                    #>gi|23455713|ref|NC_004301.1| Enterobacteria phage FI, complete genome >gi|15183|emb|X07489.1| Bacteriophage SP genomic RNA
                    headers = line.split('>')
                    gi2 = -1
                    if len(headers) > 2:
                        hdr2 = headers[2]
                        if hdr2.startswith('gi|'):
                            gi2 = int(hdr2.split('|',2)[1])
                            if not giSeen.has_key(gi2):
                                giSeen[gi2] = 0
                    (gifld,gi,accfld,acc,txt) = line[1:].split('|',4)
                    assert gifld == 'gi'
                    gi = int(gi)
                    valsDump = inpDump.readline().rstrip('\n').split(' ') #empty string fields will be ok
                    giDump = int(valsDump[0])
                    taxidDump = int(valsDump[1])
                    try:
                        lineage = self.taxaTree.getNode(taxidDump).lineageRanksStr()
                    except KeyError:
                        lineage = 'NULL'
                        print "Warning: Lineage not found for taxid %s" % (taxidDump,)
                    if giDump == gi:
                        line = ">gi|%s|%s|%s|" % (gi,accfld,acc) + \
                            ''.join(["%s:%s " % (fld,val) for (fld,val) in zip(fldsDump[1:],valsDump[1:])]) + \
                            "lineage:%s " % (lineage,) + \
                            txt
                        if gi2 > 0:
                            giSeen[gi2] |= FGI_WRITE
                        mismatchRun = 0
                    elif giDump == gi2:
                        giSeen[giDump] |= FGI_SKIP
                        skip = True
                        mismatchRun = 0
                    else:
                        if mismatchRun >= 10:
                            raise ValueError("Mismatch between FASTA and SQl DUMP input streams:" + \
                                "fastaTitle = %s valsDump = %s" % (line,' '.join(valsDump)))
                        else:
                            if not giSeen.has_key(giDump):
                                giSeen[giDump] = 0
                            giSeen[giDump] |= (FGI_SKIP | FGI_MISM)
                            skip = True
                            mismatchRun += 1
                            print "GI mismatch for ", giDump
                    if iRec % 10000 == 0:
                        print "Done %s records" % (iRec,)
                    iRec += 1
                elif skipSeq:
                    skip = True
                if not skip:
                    outFasta.write(line)
            except:
                print "Exception with input line: ", line
                pipe.close()
                raise
        pipe.close()
        inpDump.close()
        outFasta.close()
        print 'giSeen = \n', sorted(giSeen.items())
        print 'len(giSeen) = ',len(giSeen)
        for (gi,val) in sorted(giSeen.items()):
            if val & FGI_SKIP and not val & FGI_WRITE:
                print "%s never written %s" % (gi,val)
            if val & FGI_MISM:
                print "%s mistamtch %s" % (gi,val)


def fastaToHdf(fileFastaGz,fileHdf):
    hdf = tables.openFile(fileHdf, mode='w')
    a = tables.StringAtom(itemsize=1)
    # Use ``a`` as the object type for the enlargeable array.
    seq = hdf.createEArray(hdf.root,
        'seq',
        a,
        (0,),
        "Chars",
        expectedrows=7*10**9,
        filters=tables.Filters(complevel=1, complib='zlib',shuffle=False) #lzo
        )
    #size of internal memory buffer of this leaf
    seq.nrowsinbuf = 1000000
    inp = FastaReader(openGzip(fileFastaGz,'r'))
    iRec = 0
    for rec in inp.records():
        line = ''.join((l for l in rec.seqLines()))
        #for line in rec.seqLines():
        seq.append(numpy.fromstring(line,dtype='S1'))
        if iRec % 100000 == 0:
            print "Done %d FASTA records for %s characters" % (iRec,seq.nrows)
            #break
        iRec += 1
    inp.close()
    print "Converted %d characters" % (seq.nrows)
    print seq[seq.nrows-10:seq.nrows-1]
    hdf.close()
    
def hdfToFasta(fileHdf):
    hdf = tables.openFile(fileHdf, mode='r')
    seq = hdf.root.seq
    seq.nrowsinbuf = 1000000
    chunkLen = 1000000
    for iRec in xrange(seq.nrows/chunkLen):
        chunk = seq[iRec*chunkLen:(iRec+1)*chunkLen].copy()
    print "Read %d characters" % (seq.nrows)
    hdf.close()