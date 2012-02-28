#!/bin/bash
set -ex
asmDir=/usr/local/depot/projects/GOS/baltic/assembly/CombinedAssembly_b123Viral1_96ID
python $MGT_HOME/bin/454_contig_read_cnt.py < $asmDir/454ReadStatus.txt > weights.csv
mgt-icm-classifier --mode predict --inp-seq $asmDir/454LargeContigs.fna --inp-seq-attrib weights.csv --pred-min-len-samp 1000 --pred-out-dir results --run-mode batchDep --lrm-user-options "-P 0413"
#-P 9223

