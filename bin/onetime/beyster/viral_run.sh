#!/bin/bash
set -ex
asmDir=/usr/local/projects/GOSII/lzeigler/baltic/assembly/contigs
#python $MGT_HOME/bin/454_contig_read_cnt.py < $asmDir/454ReadStatus.txt > weights.csv
mgt-icm-classifier --mode predict --inp-seq $asmDir/baltic_vir_contigs.fna --pred-mode host --pred-out-dir results --run-mode batchDep --lrm-user-options "-P 0413"
#-P 9223

