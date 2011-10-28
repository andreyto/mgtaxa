#!/bin/bash
. /home/mgtaxadev/work/packages/mgtaxa/etc/mgtaxa.shrc
python $MGT_HOME/bin/454_contig_read_cnt.py < /usr/local/depot/projects/GOS/baltic/assembly/GS667_GCVIDIU02_0p1um/454ReadStatus.txt > weights.csv
mgt-icm-classifier --mode predict --inp-seq /usr/local/depot/projects/GOS/baltic/assembly/GS667_GCVIDIU02_0p1um/454LargeContigs.fna --inp-seq-attrib weights.csv --pred-min-len-samp 1000 --pred-out-dir my_results --run-mode batchDep --lrm-user-options "-P 0116"

