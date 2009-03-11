"""Application interface to graphical representation of real valued features."""

from MGT.Svm import *
from MGT.App import *
import mdp

__all__ = ["FeatGraphApp"]

def rowsDot(x):
    return (x*x).sum(1)

def coordXyzToSpherical(xyz):
    r = n.sqrt(rowsDot(xyz))
    xy = xyz[:,:-1]
    z = xyz[:,-1]
    rpt = n.zeros_like(xyz)
    rpt[:,0] = r
    S = n.sqrt(rowsDot(xy))
    #rpt[:,1] = z/r
    rpt[:,1] = n.arccos(z/r)
    x = xyz[:,0]
    y = xyz[:,1]
    #rpt[:,2] = y/S
    aY = n.arcsin(y/S)
    rpt[:,2] = n.select([x >= 0, x < 0], [aY, n.pi-aY])
    return rpt

class FeatGraphApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=("default",),
            dest="mode",default="default"),
            make_option("-p", "--method",
            action="store", type="choice",choices=("xyz",),
            dest="method",default="xyz"),
            make_option("-i", "--in-feat",
            action="append", type="string",dest="inFeat"),
            make_option("-a", "--labels",
            action="append", type="string",dest="labels"),
            make_option("-o", "--out-feat",
            action="store", type="string",dest="outFeat"),
        ]
        return Struct(usage = "Reduce dimentionality of real valued features\n"+\
                "%prog [options]",option_list=option_list)

    def doWork(self,**kw):
        opt = self.opt
        print "App options:\n", opt
        if opt.method == "xyz":
            return self.doXyz(**kw)

    def doXyz(self,**kw):
        opt = self.opt
        self.loadLabels()
        data = self.loadAsDense()
        x = data["feature"]
        if options.debug > 0:
            print "Avg feature values: %s" % (x.mean(0),)
            print "Std feature values: %s" % (x.std(0),)
        #TMP:
        x = stdScaleAndCenter(x)
        if False:
            pos = (0,1,2)
            #pos = (0,2,1) separates coastal from ocean somehow
            xn = n.zeros_like(x)
            xn[:,pos[0]] = x[:,0]
            xn[:,pos[1]] = x[:,1]
            xn[:,pos[2]] = x[:,2]
            x = coordXyzToSpherical(xn)
            x = stdScaleAndCenter(x)
        x *= (len(x)**(1./3))*3.5/x.std(0)[0] # xyz coords should be in Angstroms with VDW-like density, we scale them
        if options.debug > 0:
            print "Avg feature values after scaling: %s" % (x.mean(0),)
            print "Std feature values after scaling: %s" % (x.std(0),)
        #TMP:
        #idRec = self.idLab.getIdToRec()
        #for rec in data:
        #    rec["label"] = idRec[rec["id"]]["split"]
        data = makeDenseFeature(data["label"],x,data["id"])
        self.saveAsXyz(data)

    def loadLabels(self):
        self.idLab = loadIdLabelsMany(fileNames=self.opt.labels)
        
    def loadAsDense(self):
        sparse = loadSparseSeqsMany(self.opt.inFeat,idLab=self.idLab)
        return sparseToDenseSeqs(sparse)

    def saveAsXyz(self,x):
        saveDenseSeqsAsXyz(data=x,out=self.opt.outFeat)

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(FeatGraphApp)

