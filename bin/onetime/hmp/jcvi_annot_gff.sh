#!/bin/bash
set -ex
#for inpFasta in AssemblyRun/HMPMDA0100_allCA_withUNIQdeg.fasta; do
for inpFasta in AssemblyRun/*.fasta; do
    proj=$(basename $inpFasta)
    proj=${proj%%.fasta}
    projArea=AssemblyRun/$proj
    rm -rf annot/$proj
    python ~/work/mgtaxa/bin/jcvi_annot_gff.py $projArea/prok-annotation/camera_annotation_parser.raw.combined.out.gz $projArea/prok-annotation/input/metagene_mapped_pep.fasta.gz $projArea.fasta None 1000 annot/$proj
done
