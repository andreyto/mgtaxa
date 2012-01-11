#!/bin/bash
#rsync -av /usr/local/archive/projects/GOS/mgtaxa/db-cache/refseq-seqdb /usr/local/archive/projects/GOS/mgtaxa/db-cache/taxonomy .
#python  $MGT_HOME/MGT/ImmClassifierApp.py --mode bench  --db-seq refseq-seqdb --bench-frag-len-list 400,10000 --db-bench-frag-count-max 100 --run-mode batchDep --lrm-user-options '-P 0413'
python  -m pdb $MGT_HOME/MGT/ImmClassifierApp.py --mode proc-bench-scores  --bench-n-lev-test-model-min 3 --cwd=`pwd` --run-mode inproc

