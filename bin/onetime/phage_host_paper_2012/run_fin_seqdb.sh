#!/bin/bash
python $MGT_HOME/MGT/ImmClassifierApp.py --mode fin-ref-seqdb \
--inp-ncbi-seq $MGT_NCBI_DB'/refseq/*.fna.gz' \
--db-seq $MGT_DATA/refseq-seqdb --run-mode batchDep \
--lrm-user-options '-P 9223'

