from subprocess import check_call
from os.path import join as pjoin
import os

annot_cons_script=os.path.expandvars("$MGT_HOME/MGT/JCVI/AnnotConsensusApp.py")
annot_cons_cmd=["python",annot_cons_script,"classify"]

apis_top="/usr/local/projects/GOS3/Tier1/jbadger/omz"
orf_top="/usr/local/projects/GOS3/Tier1/annotation/omz"
for grp_apis,grp_orf in (("ctg","ctg"),("single","singleton")):
    apis_annot=pjoin(apis_top,grp_apis,"OMZ{0}_apis.json.bz2".format(grp_apis))
    orf_fasta=pjoin(orf_top,"OMZ.{0}".format(grp_orf),"orf-calling/metagene_mapped_pep.fasta")
    asm_fasta="/usr/local/archive/projects/GOS/OMZ/Assembly/CA_8/9-terminator/OMZ.{0}.fasta".format(grp_orf)
    grp_work_dir=os.path.abspath(grp_orf)
    os.mkdir(grp_work_dir)
    curdir=os.getcwd()
    try:
        os.chdir(grp_work_dir)
        check_call(annot_cons_cmd+\
                [
                    "--apisAnnot",apis_annot,
                    "--inpPep",orf_fasta,
                    "--inpContigs",asm_fasta,
                ]
                )
    finally:
        os.chdir(curdir)

