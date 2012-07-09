#!/bin/bash
set -ex
projTop=$(pwd)
proj=asm_combined_454_large_5k_annotation
projAnnot=$projTop/$proj
mkdir -p annot
rm -rf annot/$proj
python -m pdb ~/work/mgtaxa/bin/jcvi_annot_gff.py $projAnnot/camera_annotation_parser.raw.combined.out.gz None $projTop/asm_combined_454_large.5K.fna None 1000 annot/$proj

