#!/bin/bash
set -ex
#python $MGT_HOME/MGT/ImmClassifierApp.py --mode make-ref-seqdb \
mgt-icm-classifier --mode make-ref-seqdb \
--inp-ncbi-seq $MGT_NCBI_DB'/refseq/*/*.fna.gz' \
--db-seq $MGT_DATA/refseq-seqdb \
--run-mode batchDep \
--lrm-user-options '-P 9223' \
--batch-backend makeflow \
--workflow-file make-ref-seqdb.mkf

MAKEFLOW_BATCH_QUEUE_TYPE=sge BATCH_OPTIONS='-P 9223 -b n -l fast' $MGT_HOME/bin/mgt_wrapper makeflow -J 1000 make-ref-seqdb.mkf

