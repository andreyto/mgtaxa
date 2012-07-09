#!/bin/bash
set -ex
asmDir=/usr/local/depot/projects/GOS/baltic/assembly/final_combined_assembly.20120520
echo "contig00111" > contig_filt.csv
#python $MGT_HOME/bin/onetime/beyster/454_contig_lib_cnt.py $asmDir contig_filt.csv
python $MGT_HOME/bin/onetime/beyster/454_contig_lib_cnt.py $asmDir

