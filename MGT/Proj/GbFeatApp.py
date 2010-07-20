### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Extracting various features from GenBank files.
"""


from MGT.Taxa import *
from MGT.FastaIO import FastaReader
from MGT.App import *

from MGT.DirStore import *

from MGT.BioUtil import *
from MGT.Kmers import *

from glob import iglob


class GbFeatApp(App):
    """App-derived class for GenBank file manipulation"""

    #batchDepModes = ("parscan",)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        pass
   
    def getDbSql(self):
        """Allocate (if necessary) and return a connection to SQL server"""
        if not hasattr(self,"dbSql"):
            self.dbSql = DbSqlMy(db="gbfeat")
        return self.dbSql

    def delDbSql(self):
        """Call this to free a connection to SQL server if it will not be needed for extended period of time"""
        if hasattr(self,"dbSql"):
            self.dbSql.close()
            del self.dbSql

    def init(self):
        self.cvTreeExe = "/home/atovtchi/work/distros/cvtree/cvtree/cvtree"
        self.cvTreeMatExe = "/home/atovtchi/work/distros/cvtree/cvtree/batch_dist.pl"
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
    
    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "split-cds":
            return self.splitAllGenBankCds()
        if opt.mode == "make-feat":
            return self.makeFeat()
        else:
            raise ValueError("Unknown mode: %s" % (opt.mode,))

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree


    def splitAllGenBankCds(self):
        """Run splitGenBankCds for all files of a given type."""
        opt = self.opt
        store = self.store
        for orgType in opt.orgTypes:
            inpGb = openCompressed(pjoin(options.refSeqDataDir,"%s.genomic.gbff.gz" % orgType),"r")
            outCdsFile = store.getFilePath("%s.cds.fasta.gz" % orgType)
            outNcdsFile = store.getFilePath("%s.ncds.fasta.gz" % orgType)
            outCds = openCompressed(outCdsFile,"w")
            outNonCds = openCompressed(outNcdsFile,"w")
            extractCdsAndNonCds(inpGb=inpGb,outCds=outCds,outNonCds=outNonCds)
            inpGb.close()
            outCds.close()
            outNonCds.close()

    def makeFeat(self):
        """Create feature vectors out of FASTA files."""
        maxSampLen = 100000
        kmerCnt = KmerSparseFeatures(sampLen=maxSampLen,
                kmerLen=2,
                rcPolicy=RC_POLICY.MERGE,
                normPolicy=NORM_POLICY.FREQ)
        for fastaFile in self.store.getFilePaths("*.fasta.gz"):
            inpFasta = FastaReader(fastaFile)
            iRec = 0
            for rec in inpFasta.records():
                kmerCnt.process(rec.sequence(format="array"))
                iRec += 1
                if iRec > 100:
                    break
            inpFasta.close()
            feat = kmerCnt.kmerFrequencies()
            id = stripSfx(os.path.basename(fastaFile),".fasta.gz")
            print id, feat


def run_GbFeat():
    opt = Struct()
    opt.runMode = "inproc" #"batchDep"
    opt.orgTypes = ["microbial"]
    modes = ["make-feat"] #["split-cds"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = GbFeatApp(opt=opt)
        jobs = app.run(depend=jobs)

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(GbFeatApp)
