"""Driver script to find CRISPR in Sloan-Air 454 reads"""
from MGT.Common import *
from MGT.Proj.CrisprApp import *

inpFastaFiles = [
"/usr/local/projects/BIOINFO/TEMP/yellowstone/yellow1/all.fasta",
"/usr/local/projects/BIOINFO/TEMP/yellowstone/yellow2/all.fasta",
"/usr/local/projects/BIOINFO/TEMP/yellowstone/yellow3/all.fasta",
"/usr/local/projects/BIOINFO/TEMP/yellowstone/yellow4/all.fasta",
"/usr/local/projects/BIOINFO/TEMP/yellowstone/yellow5/all.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_01_scf/9-terminator/ystoneSite01.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_01_scf/9-terminator/ystoneSite01.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_01_scf/9-terminator/ystoneSite01.singleton.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_02_scf_sff/9-terminator/ystoneSite02.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_02_scf_sff/9-terminator/ystoneSite02.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_02_scf_sff/9-terminator/ystoneSite02.singleton.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_03_scf/9-terminator/ystoneSite03.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_03_scf/9-terminator/ystoneSite03.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_03_scf/9-terminator/ystoneSite03.singleton.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_04_scf_sff/9-terminator/ystoneSite04.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_04_scf_sff/9-terminator/ystoneSite04.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_04_scf_sff/9-terminator/ystoneSite04.singleton.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_08_scf/9-terminator/ystoneSite08.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_08_scf/9-terminator/ystoneSite08.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_08_scf/9-terminator/ystoneSite08.singleton.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_18_scf/9-terminator/ystoneSite18.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_18_scf/9-terminator/ystoneSite18.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_18_scf/9-terminator/ystoneSite18.singleton.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_19_scf/9-terminator/ystoneSite19.deg.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_19_scf/9-terminator/ystoneSite19.scf.fasta",
"/usr/local/archive/bioinformatics/drusch/yellowstoneCSP/DougYNP/YNP_Site_19_scf/9-terminator/ystoneSite19.singleton.fasta"
]

run_name="doug_ystone"
optB = Struct()
optB.runMode = "inproc" #"batchDep"
optB.sqlEngine = "sqlite"
optB.sqlDb = run_name
optB.sqlHost = "mgtaxa-app.jcvi.org"
topRunName = pjoin(os.environ["GOSII_WORK"],run_name)

allJobs = []

for (iFile,inpFastaFile) in enumerate(inpFastaFiles):
    inpId = "%04d-%s" % (iFile+1,stripSfx(os.path.basename(inpFastaFile),".fasta"))
    #need to make a fresh copy every loop otherwise all derived opt members will
    #get cached between loops
    opt = copy(optB)
    opt.topWorkDir = pjoin(topRunName,inpId)
    makedir(opt.topWorkDir)
    opt.cwd = opt.topWorkDir
    opt.inpFastaFiles = [ inpFastaFile ]
    opt = CrisprApp.fillWithDefaultOptions(opt)

    jobs = []
    
    opt.mode = "import-inp-seq"
    #app = CrisprApp(opt=opt)
    #print "Running mode %s" % opt.mode
    #jobs = app.run(depend=jobs)
    #continue
    
    opt.mode = "findcr"
    app = CrisprApp(opt=opt)
    print "Running mode %s" % opt.mode
    jobs = app.run(depend=jobs)
    allJobs += jobs

#modes = [ "blastcr","loadbl","annotmic","stats" ] 
#modes = [ "annotmic" ]
#modes = [ "loadmic","blastdb","pilercr","loadcr","exportcr","blastcr","loadbl","annotmic","stats" ] 
#modes = [ "stats" ] #"annotmic" "loadmicprot" "annotmic" "pilercr" "loadcr" "exportcr" "loadmic" "loadannot" "exportaclame" "exportmgt" "blastcr" "loadbl" "blastdb"

