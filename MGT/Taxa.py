"""Combines classes for processing taxonomy trees"""

from MGT.TaxaTree import *
from MGT.TaxaTreeDb import *
from MGT.TaxaIO import *
from MGT.TaxaIODb import *

def loadTaxaTree(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile):
    return TaxaTree(NodeStorageNcbiDump(ncbiDumpFile=ncbiDumpFile,
        ncbiNamesDumpFile=ncbiNamesDumpFile))

def makeGiTaxBin(ncbiDumpFiles,outFile):
    """Create and save a pickled numpy gi->taxid index from a list of ncbi dump files.
    Typically, there are two dump files: one for nucleotide and another for protein sequences.
    This function checks that no GI is present in more than one file.
    The rsulting file can be loaded back into memory with loadGiTaxBin()."""
    dst = None
    for ncbiDumpFile in ncbiDumpFiles:
        inp = openCompressed(ncbiDumpFile,'r')
        twocol = numpy.fromfile(inp,dtype="i4",sep='\n')
        twocol.shape = (twocol.shape[0]/2,2)
        gi2taxa = fromWhereItems({'ind':twocol[:,0], 'val':twocol[:,1]})
        print gi2taxa.shape, gi2taxa.dtype
        if dst is None:
            dst = gi2taxa
        else:
            if len(dst) >= len(gi2taxa):
                src = gi2taxa
            else:
                src = dst
                dst = gi2taxa
            dst_sub = dst[:len(src)]
            assert numpy.logical_not(numpy.logical_and(dst_sub > 0,src > 0)).all()
            src_ind = numpy.where(src>0)
            dst[src_ind] = src[src_ind]
    dumpObj(dst,outFile)

def loadGiTaxBin(inFile):
    return loadObj(inFile)

