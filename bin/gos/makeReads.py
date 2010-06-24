### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from MGT.Common import *
from bin.gos import *
from bin import gosToSvm

def makeSangerReads():
    makedir(GOS_VIR_SANG_DIR)
    opt,args = gosToSvm.getProgOptions(args=[])
    for inpFasta in GOS_VIR_ANNOT_FASTA_A:
        outFasta = pjoin(GOS_VIR_SANG_DIR,os.path.basename(inpFasta))
        gosToSvm.fastaToSvm(inFileFasta=inpFasta,outName=outFasta,opt=opt)

def listReadFiles():
    return glob.glob(pjoin(GOS_VIR_454_DIR,"*.fasta.seq")) + \
            glob.glob(pjoin(GOS_VIR_SANG_DIR,"*.fasta"))

if __name__ == "__main__":
    makeSangerReads()

