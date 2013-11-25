#!/bin/sh
this_dir=$(dirname $0)
python $this_dir/convert_shibu_fasta.py '/usr/local/depot/projects/GOS/analysis/syooseph/gos3asm_2012/annotation/sequenced_genomes/fasta/*.fasta' genome.mod.fasta genome.mod.json
python $this_dir/convert_shibu_fasta.py '/usr/local/depot/projects/GOS/analysis/syooseph/gos3asm_2012/annotation/taxonomy/fasta/*.fasta' scaff.mod.fasta scaff.mod.json

