### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""GOS-specific processing environment"""
import os
import os.path.join as pj
GOS_PROJ_DIR="/usr/local/projects/GOSII/"
GOS_ANNOT_DIR=pj(GOS_PROJ_DIR,"ANNOTATION")
GOS_VIR_ANNOT_DIRS=(\
        "GSIOVIR108-I-01-2-3KB_reprocessed",
        "GSIOVIR110-I-01-3-4KB_reprocessed",
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

