"""Process per-sample assemblies for Beyster GOS"""
from MGT.ImmClassifierApp import *

topSampDir = "/usr/local/depot/projects/GOS/baltic"
topAsmDir = pjoin(topSampDir,"assembly")

sampSubDirs = \
"""GS659_GCX5D7Y01_0p1um
GS665_GCX5D7Y01_0p1um
GS667_GCVIDIU02_0p1um
GS669_GCVIDIU02_0p1um
GS677_GCVEAXJ02_0p1um
GS678_GCVEAXJ02_0p1um
GS679_GDNEDKP02_3p0um
GS679_GLDFQNX01_0p1um
GS680_GDNEDKP02_3p0um
GS681_GCXE2IL02_0p1um
GS683_GCXE2IL02_0p1um
GS685_GCZC3J301_0p1um
GS687_GCZC3J301_0p1um
GS688_GCZC3J302_0p1um
GS689_GCZC3J302_0p1um
GS695_GDQ27C301_0p1um
"""

sampSubDirs = [ l.strip() for l sampSubDirs.split("\n") if l.strip() ]

print sampSubDirs
