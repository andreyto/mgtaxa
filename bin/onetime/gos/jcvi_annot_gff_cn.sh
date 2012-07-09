#!/bin/bash
set -ex
proj=/usr/local/projects/GOS3/Tier2/CN_assembly_annotation_2012/CN_2012_92pct
python ~/work/mgtaxa/bin/jcvi_annot_gff.py \
    $proj/apis/CN_GOS4_92.tab.gz \
    None \
    None \
    $proj/input/GOS3.scf.fasta \
    None \
    10 \
    annot1
