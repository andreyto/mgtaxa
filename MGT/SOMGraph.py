from MGT.Common import *
import pylab as pl
import matplotlib as plib

def somPlotScatter(mod,colorMap,markerSize=64,alpha=0.75,figSizeScatter=(16,16),figSizePie=(8,8)):
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

    pl.figure(1, figsize=figSizeScatter)
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
    pl.figure(2, figsize=figSizePie)
    labCounts = mod.sampLabelCounts()
    #pdb.set_trace()
    pl.pie(labCounts.values(), labels=labCounts.keys(), colors=[ colorMap[l] for l in labCounts.keys() ],autopct='%1.1f%%', shadow=True)

    pl.show()

