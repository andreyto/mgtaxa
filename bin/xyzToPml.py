from MGT.Svm import *
from pymol import cmd, stored

class PropertyMapper:

    objName = "x"

    def __init__(self,xyzFile,idMapFile):
        cmd.load(xyzFile,self.objName)
        stored.ids = loadSeqsIdDef(xyzFile).tolist()
        stored.idMap = loadObj(idMapFile).getIdToRec()
        cmd.alter(self.objName,"text_type=stored.ids.pop(0)")
        cmd.alter(self.objName,"chain=1")
        cmd.alter(self.objName,"segi=numeric_type")
        cmd.alter(self.objName,"label=stored.idMap[text_type]['oldId']")
        #for lab in range(10): cmd.color(lab,"(nt. %s)" % lab)
        #cmd.spectrum("b","selection")

    def colorByProp(self,prop,palette="rainbow"):
        stored.propUniqVals = set()
        cmd.iterate(self.objName,"stored.propUniqVals.add(%s)" % prop)
        v = sorted(stored.propUniqVals)
        b = n.arange(1,len(v)+1,dtype=float)# / len(v)
        stored.b = dict(zip(v,b))
        cmd.alter(self.objName,"b=stored.b[%s]" % prop)
        cmd.spectrum("b",palette,self.objName)


## pseudoatom example - usefull to add a few points with arbitrary
## e.g. radii to mark some locations, but gets very slow for 100 poitns or more
#cmd.pseudoatom("tmp",elem='Ar', b=40, color=words[-1], pos=pos,vdw=1.7)

