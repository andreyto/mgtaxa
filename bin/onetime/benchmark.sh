#!/bin/bash
#rsync -av /usr/local/archive/projects/GOS/mgtaxa/db-cache/refseq-seqdb /usr/local/archive/projects/GOS/mgtaxa/db-cache/taxonomy .
#All but viruses, all models
python  $MGT_HOME/MGT/ImmClassifierApp.py --mode bench  --db-bench-taxids-exclude-trees 10239,12884,28384,12908 --db-bench-filter-by-models 1 --db-seq refseq-seqdb --bench-frag-len-list 100 --db-bench-frag-count-max 200 --bench-n-lev-test-model-min  1 --run-mode batchDep --lrm-user-options '-P 0413'
#Viruses, viral models
#python  $MGT_HOME/MGT/ImmClassifierApp.py --mode bench  --score-taxids-exclude-trees 131567 --db-bench-filter-by-models 1 --db-seq refseq-seqdb.vir --bench-frag-len-list 100,400,1000,5000,10000 --db-bench-frag-count-max 10 --bench-n-lev-test-model-min  1 --run-mode batchDep --lrm-user-options '-P 0413'
#python  -m pdb $MGT_HOME/MGT/ImmClassifierApp.py --mode proc-bench-scores  --bench-n-lev-test-model-min 1 --db-bench-frag-len 100 --cwd=`pwd` --run-mode inproc
#python  -m pdb $MGT_HOME/MGT/ImmClassifierApp.py --mode proc-bench-scores  --bench-n-lev-test-model-min 1 --db-bench-frag-len 400 --cwd=`pwd` --run-mode inproc
#python  -m pdb $MGT_HOME/MGT/ImmClassifierApp.py --mode proc-bench-scores  --bench-n-lev-test-model-min 1 --db-bench-frag-len 400 --cwd=`pwd` --run-mode inproc --bench-out-dir benchResults.marine --bench-proc-subtrees 2.Mooore_taxid.csv --bench-proc-subtrees 3.Moore_paper_extra_taxid.csv --bench-proc-subtrees 4.Cyanobacteria_Roseobacter_SAR11_taxid.csv
#for frLen in 100 400 1000 5000 10000; do cat immClassifierApp.taI2d8.work/benchWorkDir/$frLen/benchResults/bench.aggr.csv; done
