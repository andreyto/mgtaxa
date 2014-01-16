"""Entry point for Galaxy MGTAXA tools"""

from MGT.Common import *
from MGT.RigidCli import *
from subprocess import check_call
import shlex

_mgt_common_options = """ \
--run-mode batchDep \
--batch-backend makeflow \
--workflow-run 1 """


def train(argv_dict,argv):
    argv = rigid_cli_dict_to_argv(argv_dict)
    cmd_base = shlex.split(
            """\
            mgt-icm-classifier \
            --inp-train-seq-format generic \
            --db-seq seq.db """ + \
            _mgt_common_options
            ) + \
        argv
    
    cmd = cmd_base + \
        shlex.split("--mode make-ref-seqdb")
    
    check_call(cmd)
    
    cmd = cmd_base + \
        shlex.split("--mode train")
    
    check_call(cmd)

def predict(argv_dict,argv):
    d = argv_dict
    db_imm_default = d["--db-imm-default"]
    if db_imm_default[0].lower() == "yes":
        db_imm_default[0] = "1"
    elif db_imm_default[0].lower() == "no":
        db_imm_default[0] = "0"
    argv = rigid_cli_dict_to_argv(d)
    cmd = shlex.split(
            """\
            mgt-icm-classifier \
            --db-seq seq.db """ + \
            _mgt_common_options
            ) + \
        argv + \
        shlex.split("--mode predict")
    
    check_call(cmd)


cli_patterns = {
        
        "train" : {
            "args" : """\
                    --inp-train-seq $param \
                    --inp-train-model-descr $param \
                    --db-imm-archive $param \
                    --lrm-user-options $param \
                    --makeflow-options $param \
                    --web $param""",
            "func" : train
            },

        
        "predict" : {
            "args" : """\
                    --db-imm-default $param \
                    --db-imm-archive $param \
                    --inp-seq $param \
                    --pred-out-taxa-csv $param \
                    --pred-out-stats-csv $param \
                    --pred-out-stats-html $param \
                    --pred-out-stats-pdf $param \
                    --pred-mode $param \
                    --score-taxids-exclude-trees $param \
                    --lrm-user-options $param \
                    --makeflow-options $param \
                    --web $param""",
            "func" : predict
            }
        }

def main():
    rigid_cli_dispatch(cli_patterns=cli_patterns)

if __name__ == '__main__':
    main()

