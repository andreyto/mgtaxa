import sys
from numpy import *
from MGT.Util import *
import pdb
labToId = loadObj("labelToId.pkl")
inp = open(sys.argv[1],'r')
lb = fromiter(( int(l.split(None,1)[0]) for l in inp ),dtype=int)
inp.close()
c = bincount(lb)
print zip(labToId,c)
print [ (x[0][0],x[1]) for x in ndenumerate(c) ]
minUse = int(sys.argv[2])
cUse = c.copy()
cUse[c < minUse] = 0
cUse[c >= minUse] = minUse
out = open(sys.argv[1]+'.b','w')
for l in open(sys.argv[1],'r'):
    lb = int(l.split(None,1)[0])
    if cUse[lb] > 0:
        if float(minUse)/c[lb] > random.random_sample():
            cUse[lb] -= 1
            out.write(l)
out.close()

