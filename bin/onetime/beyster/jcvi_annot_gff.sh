#!/bin/bash
set -ex
#python -m pdb ../bin/onetime/beyster/jcvi_annot_gff.py /usr/local/depot/projects/GOS/baltic/assembly/assembly_annotation/viral_5Kcontigs/viral-annotation/uniref_blastp_btab.combined.out.gz /usr/local/depot/projects/GOS/baltic/mgtaxa/2012-02-09/viral/results.euk_round_up/pred-taxa 10000 annot/btab
#python -m pdb ../bin/onetime/beyster/jcvi_annot_gff.py /usr/local/depot/projects/GOS/baltic/assembly/assembly_annotation/viral_5Kcontigs/viral-annotation/camera_annotation_parser.raw.combined.out.gz /usr/local/depot/projects/GOS/baltic/mgtaxa/2012-02-09/viral/results.euk_round_up/pred-taxa 10000 annot/camera
projArea=/usr/local/depot/projects/GOS/baltic
python -m pdb ../bin/jcvi_annot_gff.py $projArea/assembly/assembly_annotation/combined_assembly_1Kcontigs/prok-annotation/camera_annotation_parser.raw.combined.out.gz None None $projArea/mgtaxa/2012-02-23/combined/results.10K.taxa.vir/pred-taxa 10000 annot/ann_nodes

