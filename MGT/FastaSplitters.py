"""Methods that split FASTA streams by various criteria"""
from MGT.Taxa import *

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
                    meta_group = dict(
                            id=node.id,
                            taxid=node.id,
                            name=node.name
                            )
                    meta = dict(
                            id=gi,
                            gi=gi,
                            name=hdr
                            )
                    yield dict(seq=rec,meta=meta,meta_group=meta_group)

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
            node = taxaTree.getNode(rec["meta_group"]["taxid"])
            if not hasattr(node,'seq'):
                node.seq = []
            seqRec = Struct(gi=rec["meta"]["gi"])
            if storeHeader:
                seqRec.header = rec["seq"].header().strip()
            seqLen = None
            if storeSeq:
                seqRec.seq = rec["seq"].sequence()
                seqLen = len(seqRec.seq)
            if storeSeqLen:
                if seqLen is None:
                    seqLen = rec["seq"].seqLen()
                seqRec.seqLen = seqLen
            node.seq.append(seqRec)
        inpSeq.close()
    return taxMis

def splitFastaFilesByGroupId(iterFastaWithMeta,outStore):
    """@param outStore SeqDbFasta instance"""
    idSeen = set()
    lastId = None
    out = None
    try:
        for rec in iterFastaWithMeta:
            seq = rec["seq"]
            meta_group = rec["meta_group"]
            id = meta_group["id"]
            if lastId is None or lastId != id:
                if out is not None:
                    out.close()
                out = outStore.fastaWriterUncompr(id,mode="a")
                if not id in idSeen:
                    outStore.saveMetaDataById(id,meta_group)
                    idSeen.add(id)
            out.header(seq.header())
            out.seqLines(seq.seqLines())
    finally:
        if out is not None:
            out.close()

def splitFastaFilesByTaxa(inSeqs,outStore,taxaTree=None,giToTaxa=None,filt=None):
    """@param outStore SeqDbFasta instance"""
    from MGT.FastaIO import FastaReader
    if taxaTree is None:
        taxaTree = loadTaxaTree()
    if giToTaxa is None:
        giToTaxa = loadGiTaxBin()
    taxMis = Struct()
    if filt is None:
        filt = lambda x: x
    def _multi_iter():
        for inSeq in inSeqs:
            inpSeq = FastaReader(inSeq)
            for rec in ncbiFastaRecordsWithTaxa(fastaReader=filt(inpSeq),
                    taxaTree=taxaTree,
                    giToTaxa=giToTaxa,
                    errorCounter=taxMis):
                yield rec
            inpSeq.close()
    splitFastaFilesByGroupId(iterFastaWithMeta=_multi_iter(),outStore=outStore)
    return taxMis

def splitFastaFilesByModel(inSeqs,modelsMeta,outStore,
        taxaTree=None,checkTaxa=True,filt=None):
    """
    @param modelsMeta iterator of models meta data objects
    @param outStore SeqDbFasta instance
    """
    from MGT.FastaIO import FastaReader
    if checkTaxa and taxaTree is None:
        taxaTree = loadTaxaTree()
    mis = Struct()
    if filt is None:
        filt = lambda x: x
    #creat mapping from id_seq to model metadata
    seqToModelMeta = dict()
    for meta in modelsMeta:
        for id_seq in meta["ids_seq"]:
            seqToModelMeta[id_seq] = meta
        if checkTaxa:
            assert taxaTree.getNode(meta["taxid"])
    def _multi_iter():
        for inSeq in inSeqs:
            with closing(FastaReader(inSeq)) as inpSeq:
                for seq in filt(inpSeq).records():
                    id_seq = seq.getId()
                    modelMeta = seqToModelMeta[id_seq]
                    meta_group = dict(
                            id=modelMeta["id"],
                            taxid=modelMeta["taxid"],
                            name=modelMeta["name"]
                            )
                    rec = dict(
                            seq=seq,
                            meta_group=meta_group
                            )
                    yield rec
    splitFastaFilesByGroupId(iterFastaWithMeta=_multi_iter(),outStore=outStore)


