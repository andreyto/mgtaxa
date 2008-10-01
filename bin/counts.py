import sys
from numpy import *
from MGT.Util import *
try:
    labFile = sys.argv[2]
except IndexError:
    labFile = "labelToId.pkl"
labToId = loadObj(labFile)
lb = fromiter(( int(l.split(None,1)[0]) for l in open(sys.argv[1],'r') ),dtype=int)
c = bincount(lb)
print zip(labToId,c)
print [ (x[0][0],x[1]) for x in ndenumerate(c) ]

