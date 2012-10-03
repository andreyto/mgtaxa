#!/bin/bash
set -ex
asmDir=/usr/local/projects/GOSII/lzeigler/baltic/assembly/vir/combined_balVir_newbler
python $MGT_HOME/bin/454_contig_read_cnt.py < $asmDir/454ReadStatus.txt > weights.csv
mgt-icm-classifier --mode predict --inp-seq $asmDir/454LargeContigs.fna --pred-mode host --pred-out-dir results --run-mode batchDep --lrm-user-options "-P 9223"
#-P 9223

