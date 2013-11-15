#!/bin/bash
set -ex
if true; then
#cp /usr/local/projects/GOSII/atovtchi/ref.extra.seqdb/*.fasta* $MGT_DATA/refseq-seqdb/
mgt-icm-classifier --mode train \
--db-imm $MGT_DATA/icm-refseq \
--db-seq $MGT_DATA/refseq-seqdb \
--run-mode batchDep \
--batch-backend makeflow \
--workflow-file train-icm-ref.mkf
fi

MAKEFLOW_BATCH_QUEUE_TYPE=sge \
    BATCH_OPTIONS='-P 9223 -b n -l fast' \
    $MGT_HOME/bin/mgt_wrapper \
    makeflow -J 1000 -r 3 -S 30 \
    train-icm-ref.mkf


