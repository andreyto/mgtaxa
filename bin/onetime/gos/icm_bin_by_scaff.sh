#!/bin/bash
set -ex
gen_fastas="/usr/local/depot/projects/GOS/analysis/syooseph/gos3asm_2012/annotation/sequenced_genomes/fasta/*.fasta"
asm_fastas="/usr/local/depot/projects/GOS/analysis/syooseph/gos3asm_2012/annotation/taxonomy/fasta/*.fasta"
cwd=$(pwd)
#cat $gen_fastas > gen.fasta
#cat $asm_fastas > asm.fasta
gen_fasta=$cwd/gen.fasta
asm_fasta=$cwd/asm.fasta
gen_imm=$cwd/gen.imm
asm_imm=$cwd/asm.imm
gen_tree=$cwd/gen.tree
asm_tree=$cwd/asm.tree
mgt-icm-classifier --mode train --inp-train-seq $gen_fasta --db-imm $gen_imm --taxa-tree-pkl $gen_tree --taxa-tree-ncbi-dir /usr/local/projects/GOSII/atovtchi/taxonomy --run-mode batchDep --batch-backend makeflow --workflow-file $cwd/gen.mkf
makeflow -T sge -B '-P 9223 -b n -S /bin/bash' gen.mkf

