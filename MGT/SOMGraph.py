from MGT.Common import *
#import matplotlib
#matplotlib.rcParams['font.size'] = 24
import pylab as pl
import matplotlib as plib

class SOMPlot:

    def __init__(self,mod,labColorMap,markerSize=64,alpha=0.75,figSizeScatter=12,figSizePie=8,figSizeMat=12,name="som_plot",onClick=lambda event: None):
        self.mod=mod
        self.labColorMap=labColorMap
        self.markerSize=markerSize
        self.alpha=alpha
        gridAspect=float(mod.unit.shape[1])/mod.unit.shape[0]
        self.figSizeScatter=(figSizeScatter,figSizeScatter*gridAspect)
        self.figSizeMat=(figSizeMat,figSizeMat*gridAspect)
        self.figSizePie=(figSizePie,figSizePie)
        self.iFig = 0
        self.name = name
        self.onClick = onClick

    def plot(self):
        self.iFig = 0
        self.plotScatter()
        self.plotPie()
        #self.plotMat(self.mod.umat)
        #self.plotMat(self.mod.pmat)
        #self.plotMat(self.mod.usmat)
        #pl.show()

    def nextIFig(self,iFig=0):
        iFigNext = max(self.iFig+1,iFig)
        self.iFig = iFigNext
        return iFigNext

    def figFileName(self,iFig):
        return "%s-%03i.png"%(self.name,iFig)

    def plotScatter(self,iFig=0):
        iFig = self.nextIFig(iFig)
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

        fig = pl.figure(iFig, figsize=figSizeScatter)
        pl.scatter(x,y,c=c,s=s, alpha=alpha,picker=5)

        #ticks = arange(-0.06, 0.061, 0.02)
        xt = n.arange(0,unit.shape[0]+1,10,dtype=float)
        pl.xticks(xt)
        yt = n.arange(0,unit.shape[1]+1,10,dtype=float)
        pl.yticks(yt)
        pl.xlabel(r'X', fontsize=20)
        pl.ylabel(r'Y', fontsize=20)
        pl.title('Scatter plot of mapped vectors')
        pl.grid(True)
        #pl.savefig("scat.png")
        #pl.colorbar()
        pl.savefig(self.figFileName(iFig))

        cid = fig.canvas.mpl_connect('button_press_event', self.onClick)
        pl.show()


    def plotPie(self,iFig=0):
        iFig = self.nextIFig(iFig)
        pl.figure(iFig, figsize=self.figSizePie)
        labCounts = self.mod.sampLabelCounts()
        keys = sorted(labCounts.keys())
        pl.pie([ labCounts[l] for l in keys ], labels=keys, colors=[ self.labColorMap[l] for l in keys ],
                autopct='%1.1f%%', shadow=True)
        pl.savefig(self.figFileName(iFig))

    def plotMat(self,mat,iFig=0):
        iFig = self.nextIFig(iFig)
        pl.figure(iFig,figsize=self.figSizeMat)
        mode = "prod"
        if mode == "test1":
            mat = n.zeros_like(self.mod.umat)
            for x in xrange(mat.shape[0]):
                for y in xrange(mat.shape[1]):
                    mat[x,y] = x*mat.shape[1]+y*10
        else:
            pass
            #pdb.set_trace()
            #umat = umat.clip(0.03,1.)
        mat = n.transpose(mat)
        pl.imshow(mat,aspect="equal",origin="lower",interpolation="nearest")
        pl.savefig(self.figFileName(iFig))

