"""Process per-sample assemblies for Beyster GOS"""
from MGT.ImmClassifierApp import *
from MGT.Asm import *
from MGT import PostPredAn
import re

class MultiSampClassifier:

    def __init__(self):

        sampSubDirs = \
        """GS659_GCX5D7Y01_0p1um
        GS665_GCX5D7Y01_0p1um
        GS667_GCVIDIU02_0p1um
        GS669_GCVIDIU02_0p1um
        GS677_GCVEAXJ02_0p1um
        GS678_GCVEAXJ02_0p1um
        GS679_GDNEDKP02_3p0um
        GS679_GLDFQNX01_0p1um
        GS680_GDNEDKP02_3p0um
        GS681_GCXE2IL02_0p1um
        GS683_GCXE2IL02_0p1um
        GS685_GCZC3J301_0p1um
        GS687_GCZC3J301_0p1um
        GS688_GCZC3J302_0p1um
        GS689_GCZC3J302_0p1um
        GS695_GDQ27C301_0p1um
        """
        sampSubDirs = [ l.strip() for l in sampSubDirs.split("\n") if l.strip() ]
        self.sampSubDirs = sampSubDirs
        
        self.topWorkDir = os.getcwd()

        self.topPredDir = pjoin(self.topWorkDir,"predict")
        
        self.dbMergeFile = pjoin(self.topWorkDir,"pred-taxa.sqlite")
        
        self.envInpCsv = pjoin(self.topWorkDir,"baltic_meta.csv")
        
        self.topPivotDir = pjoin(self.topWorkDir,"pivot")
        
        self.csvPivotBase = pjoin(self.topPivotDir,"pred-taxa-pivot")
        
        self.topEnvDir = pjoin(self.topWorkDir,"env")
        
        self.csvEnvBase = pjoin(self.topEnvDir,"pred-taxa-env")

    def predict(self):

        topSampDir = "/usr/local/depot/projects/GOS/baltic"
        topAsmDir = pjoin(topSampDir,"assembly")


        sampDirs = [ pjoin(topAsmDir,l) for l in self.sampSubDirs ]

        #print sampSubDirs

        makedir(self.topPredDir)

        metaCsv = pjoin(self.topWorkDir,"baltic_meta.csv")

        jobsFin = []

        icmDbRef = pjoin(os.environ["GOSII_WORK"],"icm-refseq")
        for sampDir in sampDirs:
            d = os.path.basename(sampDir)
            workDir = pjoin(self.topPredDir,d)
            makedir(workDir)
            try:
                os.chdir(workDir)
                sampAttr = pjoin(workDir,"samp.attr.csv")
                outCnt = open(sampAttr,"w")
                contigReadCount454(asmDir=sampDir,out=outCnt)
                outCnt.close()

                opt = Struct()
                opt.runMode = "batchDep" #"inproc"
                opt.immDb = icmDbRef
                opt.inpSeq = pjoin(sampDir,"454AllContigs.fna")
                opt.sampAttrib = sampAttr
                opt.predMinLenSamp = 1000

                jobs = []

                for mode in ("predict",):
                    opt.mode = mode #"predict" "proc-scores" #"proc-scores-phymm" #"perf" #"proc-scores"
                    app = ImmClassifierApp(opt=opt)
                    jobs = app.run(depend=jobs)
                jobsFin += jobs
            finally:
                os.chdir(topWorkDir)
        return jobsFin

    def merge(self):
        samples = [ (sampSubDir.split("_")[0],pjoin(self.topPredDir,sampSubDir,"results/pred-taxa.sqlite")) for \
                sampSubDir in self.sampSubDirs ]
        PostPredAn.mergeStats(samples=samples,dbOut=self.dbMergeFile)
    
    def exportPivots(self):
        makedir(self.topPivotDir)
        PostPredAn.exportAsPivots(dbInp=self.dbMergeFile,csvOutBase=self.csvPivotBase,minWeight=100)

    def loadEnvData(self):
        def _hdrPreProc(row):
            ind = row.index("Latitude, Longitute (DD)")
            assert ind == 3
            return [ re.sub(r"[^a-zA-Z0-9_]+","_",name).strip("_") for name in (["id_samp"]+row[1:ind]+["Latitude (DD)", "Longitute (DD)"]+row[ind+1:]) ]
        def _preProc(row,*l,**kw):
            ind = 3
            return (row[:ind]+row[ind].split(",")+row[ind+1:],)
        PostPredAn.loadEnvData(csvInp=self.envInpCsv,dbOut=self.dbMergeFile,preProc=_preProc,hdrPreProc=_hdrPreProc)
        
    def exportEnvData(self):
        makedir(self.topEnvDir)
        PostPredAn.exportEnvData(dbInp=self.dbMergeFile,csvOutBase=self.csvEnvBase)

    def main(self):
        #self.predict()
        #self.merge()
        self.exportPivots()
        self.loadEnvData()
        self.exportEnvData()

MultiSampClassifier().main()

