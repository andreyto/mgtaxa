#!/bin/bash
#python $MGT_HOME/bin/shredFasta.py -i refseq/microbial.15.1.genomic.fna.gz -o mic_test/microbial.15.1.frag_1000.fna --frag-size 1000 --frag-count-ratio 0.05
mgt-icm-classifier --mode proc-scores \
--db-imm $MGT_DATA/icm-refseq \
--inp-seq $MGT_DATA/mic_test/microbial.15.1.frag_1000.fna \
--pred-min-len-samp 20 \
--pred-out-dir $MGT_DATA/microbial.15.1.frag_1000_results.roundup \
--run-mode inproc \
--lrm-user-options '-P 0413'

