#!/bin/bash
set -ex

PHAGE_HOST_DB=$MGT_DATA/phage-host-db

mkdir -p $PHAGE_HOST_DB
pushd $PHAGE_HOST_DB

#for mode in build-db-ph sel-db-ph-pairs; do
for mode in sel-db-ph-pairs; do

    python $MGT_HOME/MGT/PhageHostApp.py \
    --cwd `pwd` \
    --mode $mode \
    --run-mode inproc \
    --batch-backend makeflow

    #makeflow -T sge -B '-P 9223 -b n -S /bin/bash' workflow

done

