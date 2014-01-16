from MGT.Common import *
from subprocess import check_call

seqDbPath1 = pjoin(options.testDataDir,"seqdb-fasta")
seqDbPath2 = pjoin(options.testDataDir,"fasta")

lrmUserOptions = r"-b n -P 0413"
makeflowOptions = r"-T sge"

def test_train():
    cmd = """\
    python $MGT_HOME/bin/galaxy/mgtaxa.py \
    train \
    --inp-train-seq {} \
    --inp-train-model-descr {} \
    --db-imm-archive test_galaxy_tool.imm \
    --lrm-user-options "{}" \
    --makeflow-options "{}" \
    --web 0""".format(
        pjoin(seqDbPath2,"generic.mod.train.fasta.gz"),
        pjoin(seqDbPath2,"generic.mod.train.json"),
        lrmUserOptions,
        makeflowOptions
        )
    check_call(cmd,shell=True)

def test_pred():
    cmd = """\
    python $MGT_HOME/bin/galaxy/mgtaxa.py \
    predict \
    --db-imm-default 0 \
    --db-imm-archive test_galaxy_tool.imm \
    --inp-seq {} \
    --pred-out-taxa-csv test_galaxy_tool.pred.csv \
    --pred-out-stats-csv test_galaxy_tool.pred.stats.csv \
    --pred-out-stats-html test_galaxy_tool.pred.stats.html \
    --pred-out-stats-pdf test_galaxy_tool.pred.stats.pdf \
    --pred-mode taxa \
    --score-taxids-exclude-trees None \
    --lrm-user-options "{}" \
    --makeflow-options "{}" \
    --web 0""".format(
        pjoin(seqDbPath2,"generic.mod.train.fasta.gz"),
        lrmUserOptions,
        makeflowOptions
        )
    check_call(cmd,shell=True)


test_train()
test_pred()

