### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Util import openGzip, FastaReader

import tables
import numpy

import os

tmpDir = "/home/atovtchi/tmp"
#tmpDir = "/export/atovtchi"

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

fileFastaGz = "/home/atovtchi/work/phyla/phyla_sel.samp.fasta.gz"

fileHdf = os.path.join(tmpDir,"phyla_sel.samp.hdf")

#fastaToHdf(fileFastaGz,fileHdf)
hdfToFasta(fileHdf)

