#!/bin/bash
set -ex
projRoot=/usr/local/scratch/METAGENOMICS/mbihan/HMP_MDA/AssemblyRun

projBase=S_AUREUS_ToAutomon_K-MerCovRed_CLC.454Contigs
projPepBase=S_AUREUS_pep

projBase=HMPMDA0100_ToAutomon_K-MerCovRed_CLC.454Contigs
projPepBase=HMPMDA0100_pep

#projBase=E_COLIToAutomon_K-MerCovRed_Newbler.454Contigs
#projPepBase=E_COLI_pep

proj=$projRoot/$projBase
mkdir -p $projBase
cd $projBase
python ~/work/mgtaxa/bin/jcvi_annot_gff.py \
    $proj/$projPepBase.tab \
    None \
    $proj/orf-calling/$projPepBase.fasta \
    $proj/orf-calling/input/$projBase.fna \
    None \
    10 \
    annot \
    None \
    1
