#!/bin/bash
#rsync -av /usr/local/archive/projects/GOS/mgtaxa/db-cache/refseq-seqdb /usr/local/archive/projects/GOS/mgtaxa/db-cache/taxonomy .
python  $MGT_HOME/MGT/ImmClassifierApp.py --mode bench  --db-seq refseq-seqdb --bench-frag-len-list 100,400,1000,10000 --db-bench-frag-count-max 200 --bench-n-lev-test-model-min  1 --run-mode batchDep --lrm-user-options '-P 0413'
#python  -m pdb $MGT_HOME/MGT/ImmClassifierApp.py --mode proc-bench-scores  --bench-n-lev-test-model-min 1 --db-bench-frag-len 400 --cwd=`pwd` --run-mode inproc
#python  -m pdb $MGT_HOME/MGT/ImmClassifierApp.py --mode proc-bench-scores  --bench-n-lev-test-model-min 1 --db-bench-frag-len 400 --cwd=`pwd` --run-mode inproc --bench-out-dir benchResults.marine --bench-proc-subtrees 2.Mooore_taxid.csv --bench-proc-subtrees 3.Moore_paper_extra_taxid.csv --bench-proc-subtrees 4.Cyanobacteria_Roseobacter_SAR11_taxid.csv

