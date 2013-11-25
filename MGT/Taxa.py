### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Combines classes for processing taxonomy trees"""

from MGT.TaxaTree import *
from MGT.TaxaTreeDb import *
from MGT.TaxaIO import *
from MGT.TaxaIODb import *

def loadTaxaTree(ncbiDumpFile=options.taxaNodesFile,
        ncbiNamesDumpFile=options.taxaNamesFile,
        ncbiMergedDumpFile=options.taxaMergedFile,
        ncbiTaxaDumpDir=None,
        allNames=False,pklFile=None,jsonFile=None):
    if pklFile is None and jsonFile is None:
        name_trans = {
                'ncbiDumpFile':'taxaNodesFile',
                'ncbiNamesDumpFile':'taxaNamesFile',
                'ncbiMergedDumpFile':'taxaMergedFile'
                }
        if ncbiTaxaDumpDir is not None:
            dirFiles = getTaxaDumpFileNames(taxaDumpDir=ncbiTaxaDumpDir)
        else:
            dirFiles = None
        x = {}
        for name in ("ncbiDumpFile","ncbiNamesDumpFile","ncbiMergedDumpFile"):
            if dirFiles:
                x[name] = dirFiles[name_trans[name]]
            else:
                x[name] = vars()[name]
        return TaxaTree(NodeStorageNcbiDump(ncbiDumpFile=x["ncbiDumpFile"],
            ncbiNamesDumpFile=x["ncbiNamesDumpFile"],
            ncbiMergedDumpFile=x["ncbiMergedDumpFile"],
            allNames=allNames))
    elif pklFile is not None and jsonFile is None:
        storePickle = NodeStoragePickle(fileName=pklFile)
        return TaxaTree(storage=storePickle)
    elif pklFile is None and jsonFile is not None:
        storePickle = NodeStorageJson(fileName=jsonFile)
        return TaxaTree(storage=storePickle)
    else:
        raise ValueError,(pklFile,jsonFile)

def loadTaxaTreeNew(allNames=False):
    return loadTaxaTree(ncbiDumpFile=options.taxaNodesFileNew,
        ncbiNamesDumpFile=options.taxaNamesFileNew,
        ncbiMergedDumpFile=option.taxaMergedFileNew,
        allNames=allNames)

def loadTaxaTreeTest(allNames=False):
    return loadTaxaTree(ncbiDumpFile=options.taxaNodesFileTest,
        ncbiNamesDumpFile=options.taxaNamesFileTest,
        ncbiMergedDumpFile=option.taxaMergedFileTest,
        allNames=allNames)

def makeGiTaxBin(ncbiDumpFiles,outFile):
    """Create and save with numpy.save() a gi->taxid index from a list of ncbi dump files.
    Typically, there are two dump files: one for nucleotide and another for protein sequences.
    This function checks that no GI is present in more than one file.
    The resulting file can be loaded back into memory with loadGiTaxBin().
    Unused array elements are set to 0"""
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
    dumpNumpy(dst,outFile)

def loadGiTaxBin(inFile=options.taxaPickled,ncbiTaxaDumpDir=None):
    if ncbiTaxaDumpDir is not None:
        inFile = getTaxaDumpFileNames(taxaDumpDir=ncbiTaxaDumpDir)["taxaPickled"]
    return loadNumpy(inFile)

def loadGiTaxBinNew(inFile=options.taxaPickledNew):
    return loadNumpy(inFile)

