"""Process per-sample assemblies for Beyster GOS"""
from MGT.ImmClassifierApp import *
from MGT.Asm import *
from MGT import PostPredAn
import re

class MultiSampClassifier:

    def __init__(self):

        sampSubDirs = \
        """GOS669-PE-ILI5-1
        GOS677-PE-ILI6-1
        GOS731-PE-ILI4-1
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

        topSampDir = "/usr/local/projects/GOS3/Tier2/sequencing_technology_comparison/assembly_comparison"
        topAsmDir = topSampDir
        asmKmerSize=31
        asmReadLen=100


        sampDirs = [ pjoin(topAsmDir,l) for l in self.sampSubDirs ]

        #print sampSubDirs

        makedir(self.topPredDir)

        metaCsv = pjoin(self.topWorkDir,"baltic_meta.csv")

        jobsFin = []

        for sampDir in sampDirs:
            d = os.path.basename(sampDir)
            workDir = pjoin(self.topPredDir,d)
            makedir(workDir)
            try:
                os.chdir(workDir)
                asmDir = pjoin(sampDir,"velvet")
                inpFastaOrig = pjoin(asmDir,"contigs.fa.gz")
                inpFastaPred = pjoin(workDir,"pred_inp.fna")
                filterFastaByLength(inpFastaOrig,inpFastaPred,
                        minLen=300,lineLen=1000)
                
                sampAttr = pjoin(workDir,"samp.attr.csv")
                outCnt = open(sampAttr,"w")
                contigReadCountVelvet(contFasta=inpFastaPred,kmerSize=asmKmerSize,readLen=asmReadLen,out=outCnt)
                outCnt.close()
                
                opt = Struct()
                opt.runMode = "batchDep" #"inproc"
                opt.inpSeq = inpFastaPred
                opt.predMinLenSamp = 300
                opt.sampAttrib = sampAttr
                opt.predOutDir = pjoin(workDir,"results")
                opt.lrmUserOptions = '-P 9223'
                opt.mode = "predict" #"export-predictions"
                ImmClassifierApp.fillWithDefaultOptions(opt)
                jobs = []
                app = ImmClassifierApp(opt=opt)
                jobs = app.run(depend=jobs)
                #print opt

                jobsFin += jobs
            finally:
                os.chdir(self.topWorkDir)
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
        self.predict()
        #self.merge()
        #self.exportPivots()
        #self.loadEnvData()
        #self.exportEnvData()

MultiSampClassifier().main()

