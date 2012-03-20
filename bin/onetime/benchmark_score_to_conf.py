#!/usr/bin/env python
from MGT.Common import *
from MGT.ImmClassifierBench import *

benchTopPath = sys.argv[1]
dbOut = sys.argv[2]

lenSampList = [100,400,1000,10000]

dbsInp = [ (lenSamp, pjoin(benchTopPath,str(lenSamp),"benchResults/bench.sqlite")) \
        for lenSamp in lenSampList ]


benchProc = ImmClassifierBenchToScore()

benchProc.catBenchSql(dbsInp=dbsInp,dbOut=dbOut)



