### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Generate c,gamma parameters for RBF kernel grid search.
Grid is taken from libsvm practical guide"""

qsubScript = """#!/bin/tcsh
#$ -hard -P 600005
#$ -l medium -l memory=2G -l arch="lx26-eon64"
#$ -cwd
## Submit as 'qsub -b n -S /bin/tcsh -t 1-%(maxTask)i script_name'.
## Note: Apparently admins changed the default value of -b to 'y'
## and by default qstat now thinks of script_name as a binary file and does not parse it for
## embedded options (09/14/07).  Shell NEEDS to be correctly specified both at the top and (?)
## in qstat options for the user environment to be correctly sourced.
## echo "Initial environment begin"
## printenv | sort
## echo "Initial environment end"
## pstree
source $HOME/.cshrc
## printenv | sort
hostname
uname -a
pwd
top -b -n 1 | head -n 15
####
#set SVM_C_LIST=(0.05 0.08 0.1 0.12 0.2 0.5 1.0 10)
#set SVM_C_LIST=(0.0001 0.0002 0.0005 0.0008 0.001 0.002 0.005 0.01)
#set SVM_C_LIST=(0.0000001 0.000001 0.00001 0.0001)
set SVM_C_LIST=(%(cl)s)
set SVM_G_LIST=(%(gl)s)
set SVM_C=$SVM_C_LIST[$SGE_TASK_ID]
set SVM_G=$SVM_G_LIST[$SGE_TASK_ID]
set SFX="c_${SVM_C}_g_${SVM_G}"
touch *
date
source /home/atovtchi/work/mgtaxa/build/x86_64/mgtaxa.insrc.cshrc
$SVM_LIBSVM_BIN/svm-train `cat weights.log`  -c $SVM_C -g $SVM_G -e 0.01 train.svm model-$SFX.svm 
$SVM_LIBSVM_BIN/svm-predict test.svm model-$SFX.svm test-$SFX.pred 
python $MGT_HOME/bin/accuracy.py -s test.svm test-$SFX.pred > test-$SFX.perf
date
touch *
"""

from numpy import *

#ca = power(2,arange(-5,17,2,dtype='f4'))
ca = power(2,arange(-15,5,2,dtype='f4'))
ga = power(2,arange(-15,5,2,dtype='f4'))

cl = []
gl = []
for c in ca:
    for g in ga:
        cl.append(c)
        gl.append(g)

cs = ' '.join("%g" % c for c in cl)
gs = ' '.join("%g" % g for g in gl)

print qsubScript % {'cl' : cs, 'gl' : gs, 'maxTask' : len(cl)}

