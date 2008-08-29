import sys
from numpy import *
from MGT.Util import *
labToId = loadObj("labelToId.pkl")
lb = fromiter(( int(l.split(None,1)[0]) for l in open(sys.argv[1],'r') ),dtype=int)
c = bincount(lb)
print zip(labToId,c)
print [ (x[0][0],x[1]) for x in ndenumerate(c) ]

