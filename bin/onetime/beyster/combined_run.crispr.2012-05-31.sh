#!/bin/bash
set -ex
asmDir=/usr/local/depot/projects/GOS/baltic/assembly/final_combined_assembly.20120520
python -m pdb $MGT_HOME/MGT/Proj/CrisprApp.py --mode findcr --cwd crispr --inp-seq $asmDir/454LargeContigs.fna --run-mode inproc --lrm-user-options "-P 9223"

