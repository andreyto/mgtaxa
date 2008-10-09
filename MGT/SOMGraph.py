from MGT.Common import *
import pylab as pl
import matplotlib as plib

class SOMPlot:

    def __init__(self,mod,labColorMap,markerSize=64,alpha=0.75,figSizeScatter=12,figSizePie=8,figSizeUMat=12):
        self.mod=mod
        self.labColorMap=labColorMap
        self.markerSize=markerSize
        self.alpha=alpha
        gridAspect=float(mod.unit.shape[1])/mod.unit.shape[0]
        self.figSizeScatter=(figSizeScatter,figSizeScatter*gridAspect)
        self.figSizeUMat=(figSizeUMat,figSizeUMat*gridAspect)
        self.figSizePie=(figSizePie,figSizePie)

    def plot(self):
        self.plotScatter()
        self.plotPie()
        self.plotUmat()
        pl.show()

    def plotScatter(self):
        iFig=1
        mod=self.mod
        colorMap=self.labColorMap
        markerSize=self.markerSize
        alpha=self.alpha
        figSizeScatter=self.figSizeScatter
        # size in points ^2
        unit = mod.unit
        n_vect = sum( ( len(x) for x in unit.flat ) )
        x = n.zeros(n_vect,dtype=float)
        y = x.copy()
        c = n.zeros((n_vect,4),dtype=float)
        s = x.copy()
        i_vect = 0
        for (ind,vects) in n.ndenumerate(unit):
            n_v = len(vects)
            if n_v > 0:
                shifts = nrnd.uniform(low=0.2,high=0.8,size=(n_v,2))
                x[i_vect:i_vect+n_v] = shifts[:,0] + ind[0]
                y[i_vect:i_vect+n_v] = shifts[:,1] + ind[1]
                #shuffle labels because ordered would occlude one class more often
                lab = mod.sampLabels(nrnd.permutation(vects))
                col = [ colorMap[l] for l in lab ]
                c[i_vect:i_vect+n_v] = col
                s[i_vect:i_vect+n_v] = markerSize
                i_vect += n_v

        pl.figure(iFig, figsize=figSizeScatter)
        pl.scatter(x,y,c=c,s=s, alpha=alpha)

        #ticks = arange(-0.06, 0.061, 0.02)
        xt = n.arange(0,unit.shape[0]+1,dtype=float)
        pl.xticks(xt)
        yt = n.arange(0,unit.shape[1]+1,dtype=float)
        pl.yticks(yt)
        pl.xlabel(r'X', fontsize=20)
        pl.ylabel(r'Y', fontsize=20)
        pl.title('Scatter plot of mapped vectors')
        pl.grid(True)
        pl.savefig("scat.png")
        #pl.colorbar()

    def plotPie(self):
        iFig=2
        pl.figure(iFig, figsize=self.figSizePie)
        labCounts = self.mod.sampLabelCounts()
        pl.pie(labCounts.values(), labels=labCounts.keys(), colors=[ self.labColorMap[l] for l in labCounts.keys() ],
                autopct='%1.1f%%', shadow=True)

    def plotUmat(self):
        iFig=3
        pl.figure(iFig,figsize=self.figSizeUMat)
        mode = "prod"
        if mode == "test1":
            umat = n.zeros_like(self.mod.umat)
            for x in xrange(umat.shape[0]):
                for y in xrange(umat.shape[1]):
                    umat[x,y] = x*umat.shape[1]+y*10
        else:
            umat = self.mod.umat
            #pdb.set_trace()
            #umat = umat.clip(0.03,1.)
        umat = n.transpose(umat)
        pl.imshow(umat,aspect="equal",origin="lower",interpolation="nearest")

