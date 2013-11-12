#!/bin/bash
MAKEFLOW_BATCH_QUEUE_TYPE=sge BATCH_OPTIONS='-P 9223 -b n -l fast' nohup $MGT_HOME/bin/mgt_wrapper  makeflow -J 800 -r 3 -S 3 download_refseq.internal.mkf &> makeflow.log&
