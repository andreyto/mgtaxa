#!/bin/bash
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##



PRED_ROOT=~/scratch
RES_BAK=~/work/mgtdata/res

cd $PRED_ROOT

for d in pred-*; do
    RES_SRC=$d/0000001
    mkdir -p $RES_BAK/$d
    tail -n 1 $RES_SRC/*.perf > $RES_BAK/$d/all.perf
    [ -f $RES_SRC/Readme ] && cp $RES_SRC/Readme $RES_BAK/$d/
    [ -f $RES_SRC/train.qsub ] && cp $RES_SRC/train.qsub $RES_BAK/$d/
    [ -f $RES_SRC/train.ll ] && cp $RES_SRC/train.ll $RES_BAK/$d/
    [ -f $RES_SRC/train.dv ] && cp $RES_SRC/train.dv $RES_BAK/$d/
    [ -f $RES_SRC/accushell.log ] && cp $RES_SRC/accushell.log $RES_BAK/$d/
done

