"""Driver script to find CRISPR in Sloan-Air 454 reads"""
from MGT.Common import *
from MGT.Proj.CrisprApp import *

inpFastaFiles = glob.glob("/usr/local/projects/sloan/seperateEukReads/*_schmidt_nonEuk.fa")

run_name="crispr_sloan_air"
opt = Struct()
opt.runMode = "inproc" #"batchDep"
opt.sqlDb = run_name
opt.sqlHost = "mgtaxa-app.jcvi.org"
opt.topWorkDir = pjoin(os.environ["GOSII_WORK"],run_name)
opt = CrisprApp.fillWithDefaultOptions(opt)
#rmrf(opt.crArrSeqDir)
makedir(opt.crArrSeqDir)
for (iFile,inpFastaFile) in enumerate(inpFastaFiles):
    editSymlink(inpFastaFile,pjoin(opt.crArrSeqDir,"%s.fasta-%04d" % (os.path.basename(inpFastaFile),iFile)))

jobs = []
opt.mode = "exportcr"
app = CrisprApp(opt=opt)
print "Running mode %s" % opt.mode
jobs = app.run(depend=jobs)

#modes = [ "blastcr","loadbl","annotmic","stats" ] 
#modes = [ "annotmic" ]
#modes = [ "loadmic","blastdb","pilercr","loadcr","exportcr","blastcr","loadbl","annotmic","stats" ] 
#modes = [ "stats" ] #"annotmic" "loadmicprot" "annotmic" "pilercr" "loadcr" "exportcr" "loadmic" "loadannot" "exportaclame" "exportmgt" "blastcr" "loadbl" "blastdb"

