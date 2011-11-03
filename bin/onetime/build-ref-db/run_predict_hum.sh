#!/bin/bash
#python $MGT_HOME/bin/shredFasta.py -i hum_test/chr19.fna -o hum_test/chr19.frag.fna --frag-size 1000
mgt-icm-classifier --mode predict \
--db-imm $MGT_DATA/icm-refseq \
--inp-seq $MGT_DATA/hum_test/chr19.frag.fna \
--pred-min-len-samp 20 \
--pred-out-dir $MGT_DATA/hum_test.results.chr19.fr_1000.round-up \
--run-mode batchDep \
--lrm-user-options '-P 0413'

