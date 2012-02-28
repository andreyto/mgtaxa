#!/bin/bash
set -ex

if false; then
mgt-icm-classifier --mode make-ref-seqdb \
   --inp-ncbi-seq 'refseq/viral.*genomic.fna.gz' \
   --train-min-len-samp 20000 \
   --db-seq refseq-seqdb.vir --run-mode inproc \
   --lrm-user-options '-P 0413'
else
mgt-icm-classifier --mode train \
   --train-min-len-samp 2000 \
   --train-min-len-samp-model 20000 \
   --db-seq refseq-seqdb.vir --run-mode batchDep \
   --db-imm icm-refseq.vir \
   --lrm-user-options '-P 0413'
fi
