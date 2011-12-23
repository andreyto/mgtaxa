"""Benchmark support for ImmClassifierApp"""

from MGT.DirStore import *
from MGT.FastaIO import *

class ImmClassifierBenchmark(DirStore):

    fastaSfx = ".fna"
    
    def __init__(self,*l,**kw):
        DirStore.__init__(self,*l,**kw)
        self.seqDb = None
    
    def getFastaPath(self,idDb):
        return self.getFilePath("%s%s" % (idDb,self.fastaSfx))

    def listDbIds(self,iterPaths=None):
        """List DB IDs either from this store or from the externally provided iterable"""
        return list(self.fileNames(pattern="*"+self.fastaSfx,sfxStrip=self.fastaSfx,iterPaths=iterPaths))
    
    def selectIdsDb(self,immDbs):
        """Pick those IDs from SeqDb that will be used to build the benchmark.
        @param immDbs sequence of ImmDb instances - this is used to pick only those entries
        in SeqDb that also have models built against them. The self-models will
        still be excluded during testing - using ImmDb is just an easy way to
        filter out from SeqDb viruses as it was done for ImmDb"""
        immIds = set()
        for immDb in immDbs:
            immIds |= set((str(x) for x in immDb.listImmIds()))
        seqIds = set((str(x) for x in self.seqDb.getIdList()))
        return list(immIds & seqIds)

    def makeSample(self,idDb,fragSize,fragCountMax):
        """Make a sample FASTA file for a given SeqDb ID"""
        outFasta = self.getFastaPath(idDb)
        self.shredFasta(idDb=idDb,outFasta=outFasta,fragSize=fragSize,
                fragCountMax=fragCountMax)

    def catSamples(self,outFile,idsDb=None):
        """Concatenate all sample files for a list of IDs into a single file.
        @note This cats entire files, and thus relies on the newline symbol
        being present at the end of every file. The FastaIO module that is
        used to create individual files always terminates with a newline. """
        import shutil
        
        if idsDb is None:
            idsDb = self.listDbIds()
        if not hasattr(outFile,"write"):
            outFile = openCompressed(outFile,"w")
            outClose = True
        else:
            outClose = False
        for idDb in idsDb:
            inpFile = openCompressed(self.getFastaPath(idDb),"r")
            shutil.copyfileobj(inpFile,outFile,1024*1024)
            inpFile.close()
        if outClose:
            outFile.close()

    def shredFasta(self,idDb,outFasta,fragSize,
            fragCountMax,lineLen=80,outMode="w"):
        """Shred each record in multi-FASTA file into multiple records of fixed size"""
        from MGT.FeatIO import LoadSeqPreprocShred
        seqDb = self.seqDb
        seqLenTot = seqDb.seqLengths(idDb)["len"].sum()
        seqLenRatio = (fragCountMax*fragSize)/float(seqLenTot)

        if seqLenRatio < 1.:
            #PreprocShredder treats sampNum()=0 as "all samples", so we need to set it
            #to at least 1, otherwise for catLen in a huge genome we will be getting
            #all samples. Maybe shredder's behaviour should be changed to sampNum<0 => all.
            #@todo Pick coords on a virtual concatenation of all sequences, otherwise
            #it will never be quite right.
            sampNum = lambda lab,seq,id: int(rndRound(len(seq)*seqLenRatio/fragSize))
        else:
            sampNum = 0
        inpSeq = seqDb.fastaReader(idDb)
        outSeq = FastaWriter(out=outFasta,lineLen=lineLen,mode=outMode)
        shredder = LoadSeqPreprocShred(sampLen=fragSize,
                sampNum=sampNum,
                makeUniqueId=False,
                sortByStarts=True)
        catLen = min(max(seqLenTot/(fragCountMax/10),fragSize*fragCountMax*10),10**8)
        seqCat = []
        seqCatLen = 0
        seqCatStart = 0
        indFrag = 0

        def _out_frags(seqCat,ind,seqCatStart,
                shredder=shredder,idDb=idDb,outSeq=outSeq):
            seqCat = "N".join(seqCat)
            labFr,seqFr,idFr = shredder(0,seqCat,idDb)
            startsFr = shredder.getLastSampStarts()
            for iF,sF in enumerate(seqFr):
                stF = startsFr[iF]
                idF = "%s_%s pos=%s nincat=%s len=%s" % (idDb,ind+iF,seqCatStart+stF,iF,len(sF))
                outSeq.record(header=idF,sequence=sF)
            return ind+len(seqFr)
        
        for rec in inpSeq.records():
            hdr = rec.header()
            id = rec.getId()
            seq = rec.sequence()
            seqCat.append(seq)
            seqCatLen += len(seq)
            if seqCatLen >= catLen:
                indFrag = _out_frags(seqCat,indFrag,seqCatStart)
                seqCat = []
                seqCatStart += seqCatLen
                seqCatLen = 0

        if seqCatLen > 0:
            indFrag = _out_frags(seqCat,indFrag,seqCatStart)

        inpSeq.close()
        outSeq.close()

