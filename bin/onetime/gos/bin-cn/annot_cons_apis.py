#!/usr/bin/env python

from subprocess import check_call
from os.path import join as pjoin
import os

annot_cons_script=os.path.expandvars("$MGT_HOME/MGT/JCVI/AnnotConsensusApp.py")

annot_cons_cmd=["python",annot_cons_script,"classify"]
annot_sort_apis_cmd=["python",annot_cons_script,"sort-apis"]

apis_annot="/usr/local/projects/GOSII/jbadger/gos_cn_assembly/metagene_seq_apis_reclassify.json.bz2"

apis_taxa_tree=os.path.abspath("usedTaxa.1.075_no_orphan.json")

asm_fasta="/usr/local/projects/GOS3/Tier2/CN_assembly_annotation_2012/CN_2012_92pct/input/GOS3.scf.fasta"

inp_pep="/usr/local/projects/GOS3/Tier2/CN_assembly_annotation_2012/CN_2012_92pct/orf-calling/metagene_seq.faa"

apis_annot_sorted=os.path.abspath("sorted_"+os.path.basename(apis_annot))
#apis_annot_sorted=os.path.abspath("test_apis.json")

if False:
    check_call(annot_sort_apis_cmd+\
            [
                "--inpPep",inp_pep,
                apis_annot,
                apis_annot_sorted
            ]
            )
if True:
    check_call(annot_cons_cmd+\
            [
                "--apisAnnot",apis_annot_sorted,
                "--inpContigs",asm_fasta,
                "--taxaTreeJson",apis_taxa_tree
            ]
            )

