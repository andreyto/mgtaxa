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
        allNames=False):
    return TaxaTree(NodeStorageNcbiDump(ncbiDumpFile=ncbiDumpFile,
        ncbiNamesDumpFile=ncbiNamesDumpFile,
        allNames=allNames))

def loadTaxaTreeNew(allNames=False):
    return loadTaxaTree(ncbiDumpFile=options.taxaNodesFileNew,
        ncbiNamesDumpFile=options.taxaNamesFileNew,
        allNames=allNames)

def makeGiTaxBin(ncbiDumpFiles,outFile):
    """Create and save a pickled numpy gi->taxid index from a list of ncbi dump files.
    Typically, there are two dump files: one for nucleotide and another for protein sequences.
    This function checks that no GI is present in more than one file.
    The resulting file can be loaded back into memory with loadGiTaxBin()."""
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

def loadGiTaxBin(inFile=options.taxaPickled):
    return loadObj(inFile)

def loadGiTaxBinNew(inFile=options.taxaPickledNew):
    return loadObj(inFile)

def ncbiFastaRecordsWithTaxa(fastaReader,taxaTree,giToTaxa,errorCounter):
    errorCounter.bounds=0
    errorCounter.zeroG=0
    errorCounter.trN=0
    errorCounter.trV=0
    for rec in fastaReader.records():
        hdr = rec.header()
        gi = rec.getNCBI_GI()
        if len(giToTaxa) <= gi:
            errorCounter.bounds += 1
            print "giToTaxa bounds: "+hdr
        else:
            taxid = giToTaxa[gi]
            if taxid == 0:
                errorCounter.zeroG += 1
                print "zero giToTaxa: "+hdr
            else:
                try:
                    node = taxaTree.getNode(taxid)
                except KeyError:
                    errorCounter.trN += 1
                    print "no node %s %s" % (taxid,hdr)
                else:
                    yield Struct(seq=rec,node=node,gi=gi)

def mapFastaRecordsToTaxaTree(inSeqs,taxaTree,giToTaxa,
        storeHeader=False,storeSeq=False,storeSeqLen=False):
    from MGT.FastaIO import FastaReader
    if taxaTree is None:
        taxaTree = loadTaxaTree()
    if giToTaxa is None:
        giToTaxa = loadGiTaxBin()
    taxMis = Struct()
    for inSeq in inSeqs:
        inpSeq = FastaReader(inSeq)
        for rec in ncbiFastaRecordsWithTaxa(fastaReader=inpSeq,
                taxaTree=taxaTree,
                giToTaxa=giToTaxa,
                errorCounter=taxMis):
            node = rec.node
            if not hasattr(node,'seq'):
                node.seq = []
            seqRec = Struct(gi=rec.gi)
            if storeHeader:
                seqRec.header = rec.seq.header().strip()
            seqLen = None
            if storeSeq:
                seqRec.seq = rec.seq.sequence()
                seqLen = len(seqRec.seq)
            if storeSeqLen:
                if seqLen is None:
                    seqLen = rec.seq.seqLen()
                seqRec.seqLen = seqLen
            node.seq.append(seqRec)
        inpSeq.close()
    return taxMis

def splitFastaFilesByTaxa(inSeqs,taxaTree,giToTaxa,outDir):
    from MGT.FastaIO import FastaReader
    if taxaTree is None:
        taxaTree = loadTaxaTree()
    if giToTaxa is None:
        giToTaxa = loadGiTaxBin()
    taxMis = Struct()
    lastTaxid = None
    out = None
    for inSeq in inSeqs:
        inpSeq = FastaReader(inSeq)
        for rec in ncbiFastaRecordsWithTaxa(fastaReader=inpSeq,
                taxaTree=taxaTree,
                giToTaxa=giToTaxa,
                errorCounter=taxMis):
            node = rec.node
            seq = rec.seq
            taxid = node.id
            if lastTaxid is None or lastTaxid != taxid:
                if out is not None:
                    out.close()
                out = openCompressed(pjoin(outDir,"%s.fasta.gz" % taxid),"a")
            out.write(seq.header())
            for line in seq.seqLines():
                out.write(line)
        inpSeq.close()
    if out is not None:
        out.close()
    return taxMis

