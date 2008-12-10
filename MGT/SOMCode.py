"""Interface to SOMCode - an OpenSource library for building Self-Organizing Map simulators.
http://somcode.cpatc.embrapa.br
"""

from MGT.Common import *
from MGT.SOM import *
import pdb

class SOMCode:
    def __init__(self,name):
        self.name = name
        self.ivFile = name+'.somcode.iv'
        self.setModelDir(".")

    def setModelDir(self,dirName):
        makedir(dirName)
        self.modDir = dirName
    
    def getName(self,kind="unit"):
        return os.path.join(self.modDir,self.name+"."+kind)

    def writeInput(self,data):
        m = data['feature']
        nVec = m.shape[0]
        nFeat = m.shape[1]
        out = open(self.ivFile,'w')
        out.write("%s\n" % nFeat)
        frmt = "%f "*nFeat + "\n"
        for iVec in xrange(nVec):
            out.write(frmt % (tuple(m[iVec])))
        out.close()


    def writeInputFromShogunSparse(self,feat,lab):
        f = feat.get_full_feature_matrix().transpose()
        self.writeInput(data=dict(feature=f,label=lab))

    def loadWeights(self):        
        #The first line is like:
        #5 10 15  hexa gaussian
        #where 5 is input space dimensionality; 10 15 - map size
        inp = open(self.getName(kind="wgt"),'r')
        line = inp.next().strip()
        words = line.split()
        nFeat,xdim,ydim,topol,neighb = int(words[0]),int(words[1]),int(words[2]),words[3],words[4]
        assert topol == "rect","Only rectangular grid currently supported"
        grid = n.zeros((xdim*ydim,nFeat),dtype='f8')
        iLine = 0
        for line in inp:
            grid[iLine] = n.fromstring(line.strip(),dtype='f8',sep='\t')
            iLine += 1
        #grid = n.loadtxt(inp,dtype='f8')
        grid.shape = (xdim,ydim,nFeat)
        inp.close()
        return grid

    def loadModel(self,components=("weights",)):
        mod = SOMModel()
        if "weights" in components:
            weights = self.loadWeights()
            mod.weights = weights
        return mod
