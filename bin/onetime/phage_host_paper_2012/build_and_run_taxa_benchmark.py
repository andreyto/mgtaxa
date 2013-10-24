#!/bin/bash
#everything with models, no viral models
python  $MGT_HOME/MGT/ImmClassifierApp.py \
        --mode bench \
        --score-taxids-exclude-trees 10239 \
        --db-bench-filter-by-models 1 \
        --db-imm $MGT_DATA/icm-refseq \
        --db-seq $MGT_DATA/refseq-seqdb \
        --bench-frag-len-list 100,400,1000,5000,10000 \
        --db-bench-frag-count-max 100 \
        --bench-n-lev-test-model-min  1 \
        --run-mode batchDep --lrm-user-options '-P 9223' \
        --batch-backend makeflow \
        --cwd bench
