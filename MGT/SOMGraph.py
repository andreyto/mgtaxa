from MGT.Common import *
import pylab as pl

def somPlotScatter(mod):
    # size in points ^2
    unit = mod.unit
    n_vect = sum( ( len(x) for x in unit.flat ) )
    x = n.zeros(n_vect,dtype=float)
    y = x.copy()
    c = x.copy()
    s = x.copy()
    i_vect = 0
    for (ind,vects) in n.ndenumerate(unit):
        n_v = len(vects)
        shifts = nrnd.uniform(low=0.2,high=0.8,size=(n_v,2))
        x[i_vect:i_vect+n_v] = shifts[:,0] + ind[0]
        y[i_vect:i_vect+n_v] = shifts[:,1] + ind[1]
        #shuffle labels because ordered would occlude one class more often
        c[i_vect:i_vect+n_v] = nrnd.permutation(vects) 
        s[i_vect:i_vect+n_v] = 64
        i_vect += n_v

    pl.scatter(x,y,c=c,s=s, alpha=0.9)

    #ticks = arange(-0.06, 0.061, 0.02)
    xt = n.arange(0,unit.shape[0]+1,dtype=float)
    pl.xticks(xt)
    yt = n.arange(0,unit.shape[1]+1,dtype=float)
    pl.yticks(yt)

    pl.xlabel(r'X', fontsize=20)
    pl.ylabel(r'Y', fontsize=20)
    pl.title('Scatter plot of mapped vectors')
    pl.grid(True)

    pl.show()

