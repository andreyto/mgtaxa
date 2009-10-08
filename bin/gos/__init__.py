### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""GOS-specific processing environment"""

import os
from os.path import join as pj

GOS_PROJ_DIR="/usr/local/projects/GOSII"

GOS_ANNOT_DIR=pj(GOS_PROJ_DIR,"ANNOTATION")

GOS_VIR_ANNOT_DIRS=(\
        "GSIOVIR108-I-01-2-3KB_reprocessed",
        # bact contaminated as per Shannon W.: "GSIOVIR110-I-01-3-4KB_reprocessed",
        "GSIOVIR112-I-01-3-4KB",
        "GSIOVIR112-I-02-4-6KB",
        "GSIOVIR117-I-01-2-3KB",
        "GSIOVIR122-I-01-3-4KB"
        )
GOS_VIR_ANNOT_FASTA_A = [ pj(GOS_ANNOT_DIR,d,d.split('_repr')[0]+'.fasta') for d in GOS_VIR_ANNOT_DIRS ]

GOS_MIC_ANNOT_DIRS=(\
       "GSIOLG108-G-01-2-4KB",
       "GSIOLG110-G-01-3-4KB",
       "GSIOLG112-G-01-3-4KB",
       "GSIOLG117-G-01-3-4KB",
       "GSIOLG122-G-01-3-4KB",
       "GSIOLG48-G-01-8-10KB",
       "GSIOSM048-G-01-8-10KB",
       "GSIOSM108-G-01-3-4KB",
       "GSIOSM110-G-01-3-4KB",
       "GSIOSM112-G-01-3-4KB",
       "GSIOSM117-G-01-3-4KB",
       "GSIOSM122-G-01-3-4KB"
       )
GOS_MIC_ANNOT_FASTA_A = [ pj(GOS_ANNOT_DIR,d,d+'.fasta') for d in GOS_MIC_ANNOT_DIRS ]

GOS_WORK_DATA_DIR=pj(GOS_PROJ_DIR,"atovtchi")

GOS_VIR_454_DIR=pj(GOS_WORK_DATA_DIR,"reads-orig/vir/454")
# see /usr/local/projects/GOSII/atovtchi/meta/LibInfo_454.xls for correspondence with library names

# We actually hid the 110 files under "hide" subdir, so there is no need for explicit file list anymore
if False:
    _GOS_VIR_454_FASTA_BASE=(\
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FRDLTJ301.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FRDNTY201.fasta.seq",
        # bact contaminated lib 110 as per Shannon W.: "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FRHLGZN01.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FRI4PDL01.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FRJHH4T01.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FRJHH4T02.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FSP4E4R01.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FUZIQOW02.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FW651BE02.fasta.seq",
        "jtc_sff_1119619603003_Global-Sampling-Project_Viral-fraction-sequencing_FWKLEGN02.fasta.seq"
        )
    GOS_VIR_454_FASTA_A = [ pj(GOS_VIR_454_DIR,f) for f in GOS_VIR_454_FASTA_BASE ]

GOS_VIR_SANG_DIR=pj(GOS_WORK_DATA_DIR,"reads-orig/vir/sang")

