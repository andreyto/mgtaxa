"""Operations on Self Organizing Maps.
The Maps are trained by specific implementation defined elsewhere.
Here we implement some common post-processing tasks such as
U-matrix computation, segmentation and clustering on the maps."""

from MGT.Common import *

class SOMModel(Struct):

    # Pareto radius for eucl distance between length-normalized uniformly distributed samples
    paretoRadiusLenNormSamp = 0.398
    
    def setLabels(self,label):
        self.label = label

    def sampLabels(self,sampIds):
        lab = self.label
        return numpy.asarray([ lab[id] for id in sampIds ])

    def sampIds(self):
        unit = self.unit
        return n.concatenate([x for x in unit.flat])

    def sampLabelCounts(self):
        """Return dict label->count accumulated over all vectors"""
        lab = self.label
        return binCount([ lab[id] for id in self.sampIds() ])

    def makeUMatrix(self):
        w = self.weights
        # the result, accumulate sum of distances to neighbours for each cell
        u = n.zeros(w.shape[:2],dtype='f8')
        # accumulate counts of additions made to each cell in 'u'
        c = n.zeros(w.shape[:2],dtype='i4')
        #dx = n.zeros((w.shape[0]-1,w.shape[1]  ),dtype='f8')
        #dy = n.zeros((w.shape[0]  ,w.shape[1]-1),dtype='f8')
        #dd = n.zeros((w.shape[0]-1,w.shape[1]-1),dtype='f8')
        # for each cell in 'w' where it is appropriate, 
        # we store a distance to:
        # - upper neighbour
        dx = n.sqrt(((w[1:, :,:] - w[:-1,:  ,:])**2).sum(axis=2))
        # - left neighbour
        dy = n.sqrt(((w[: ,1:,:] - w[:  ,:-1,:])**2).sum(axis=2))
        # - upper left diagonal neighbour
        dd = n.sqrt(((w[1:,1:,:] - w[:-1,:-1,:])**2).sum(axis=2))
        # - upper right diagonal neighbour
        ee = n.sqrt(((w[1:, :-1,:] - w[:-1,1:,:])**2).sum(axis=2))
        # for every pair of corresponding neighbours,
        # we accumulate distance sums in 'u', and counts in 'c'
        u[  :-1,  :  ] += dx
        c[  :-1,  :  ] += 1
        u[ 1:  ,  :  ] += dx
        c[ 1:  ,  :  ] += 1
        u[  :  ,  :-1] += dy
        c[  :  ,  :-1] += 1
        u[  :  , 1:  ] += dy
        c[  :  , 1:  ] += 1
        u[  :-1,  :-1] += dd
        c[  :-1,  :-1] += 1
        u[ 1:  , 1:  ] += dd
        c[ 1:  , 1:  ] += 1
        u[  :-1, 1:  ] += ee
        c[  :-1, 1:  ] += 1
        u[ 1:  ,  :-1] += ee
        c[ 1:  ,  :-1] += 1
        # find the average in 'u'
        # @todo If learning neighbourhood was declining too fast, it seems that some 
        # nodes w/o samples assigned to
        # them might be left with weight vectors close to their initial values, leading
        # to a very noisy U-matrix
        #u[c>0] /= c[c>0]
        assert (c>0).all()
        u /= c
        #u /= (n.sqrt((w**2).sum(axis=2))/w.shape[2])
        self.umat = u

    def makeUMatrixTest(self):
        w = self.weights
        #wl = n.sqrt((w**2).sum(axis=2))
        #w /= wl[:,:,n.newaxis]
        # the result, accumulate sum of distances to neighbours for each cell
        u = n.zeros(w.shape[:2],dtype='f8')
        for i in xrange(u.shape[0]):
            for j in xrange(u.shape[1]):
                li = max(i-1,0)
                lj = max(j-1,0)
                ui = min(i+1,u.shape[0]-1)
                uj = min(j+1,u.shape[1]-1)
                d = 0.
                c = 0
                for k in xrange(li,ui+1):
                    for l in xrange(lj,uj+1):
                        if not (i == k and j == l):
                            d += n.sqrt(((w[i,j]-w[k,l])**2).sum())
                            c += 1
                d /= c
                u[i,j] = d
        self.umat = u

    def setSamples(self,samp):
        self.samp = samp

    def getSamples(self):
        return self.samp

    def mapSamples(self,pmRadius=None):
        samp = self.samp['feature']
        assert len(samp.shape) == 2 and samp.shape[1] == self.weights.shape[-1]
        if pmRadius is None:
            pmRadius = self.paretoRadius()
        nSamp = len(samp)
        w = self.weights
        gr_dim = w.shape[:-1]
        nGrid = n.multiply.accumulate(gr_dim)[-1]
        # view w as flattened except the last dimension
        wf = w.reshape(nGrid,-1)
        # coords of best node for each sample
        bn = n.zeros((nSamp,len(gr_dim)),dtype='i4')
        pmat = n.zeros(gr_dim,dtype='f8')
        pmatf = pmat.reshape(nGrid)
        pmRadiusSq = pmRadius**2
        for iSamp in xrange(nSamp):
            df = n.square(wf-samp[iSamp]).sum(axis=1)
            pmatf[df<=pmRadiusSq] += 1
            bn[iSamp] = n.unravel_index(df.argmin(),gr_dim)
            if iSamp % 100 == 0:
                print "Mapped %d samples out of %d" % (iSamp+1,nSamp)
        self.sampNode = bn
        self.pmat = pmat

    def makeUStarMatrix(self):
        pmat = self.pmat
        iord = pmat.argsort(axis=None) #sort flattened view
        plowmat = n.zeros_like(pmat)
        plowmatf = plowmat.reshape(len(iord))
        n_ord = len(iord)
        plowmatf[:] = n_ord - 1
        plowmatf[(iord,)] -= n.arange(n_ord,dtype='i4')
        plowmatf /= n_ord
        self.plowmat = plowmat
        self.usmat = self.umat*plowmat

    def _makeMat(self,dtype='f8',extraDim=tuple()):
        gr_dim = self.weights.shape[:-1]
        m = n.zeros(gr_dim+extraDim,dtype=dtype)
        nGrid = n.multiply.accumulate(gr_dim)[-1]
        mf = m.reshape((nGrid,)+extraDim)
        return Struct(m=m,mf=mf)

    def _makeUnit(self,indType):
        sampNode = self.sampNode
        gr_dim = self.weights.shape[:-1]
        grid = n.zeros(gr_dim,dtype='O')
        nGrid = n.multiply.accumulate(gr_dim)[-1]
        gridf = grid.reshape(nGrid)
        for i in xrange(nGrid):
            gridf[i] = []
        samp = self.samp
        for i in xrange(len(samp)):
            samp_n = sampNode[i]
            if indType == 'id':
                x = samp[i]['id']
            elif indType == 'ind':
                x = i
            else:
                raise ValueError("Unknown indType value %s" % indType)
            grid[samp_n[0],samp_n[1]].append(x)
        for i in xrange(nGrid):
            gridf[i] = n.asarray(gridf[i],dtype='O')
        if indType == 'id':
            self.unit = grid
        elif indType == 'ind':
            self.unitI = grid

    def makeUnit(self):
        self._makeUnit('id')
        self._makeUnit('ind')


    def makeData(self):
        self.makeUMatrix()
        self.makeUStarMatrix()

    def paretoRadius(self,maxNSamp=500,paretoRatio=0.18):
        # for k-mer frequency vectors on a unit sphere,
        # we get 0.398
        samp = self.samp['feature']
        if len(samp) > maxNSamp:
            samp = samp[(random.sample(xrange(len(samp)),maxNSamp),)]
        nSamp = len(samp)
        d = n.zeros(nSamp*(nSamp-1)/2,dtype='f4')
        c = 0
        for i in xrange(1,nSamp):
            x = samp[i]
            for j in xrange(0,i):
                d[c] = n.sqrt(n.square(x-samp[j]).sum())
                c+=1
        cnt,bin = n.histogram(d,bins=10000)
        cntsum = cnt.cumsum().astype('f4')/cnt.sum()
        paretoRad = bin[n.where(cntsum>=paretoRatio)[0][0]]
        return paretoRad

    def makeLabelMatrix(self):
        labcnt = self.sampLabelCounts()
        lab = sorted(labcnt.keys())
        labToIdLab = dict(zip(lab,range(len(lab))))
        lmat = self._makeMat(dtype='i4',extraDim=(len(labToIdLab),))
        bn = self.sampNodes
        samp = self.samp
        sampLab = self.sampLabels(samp['id'])
        sampIdLab = n.asarray([ labToIdLab[l] for l in sampLab ],dtype='i4')
        for iSamp in xrange(len(samp)):
            samp_coor = bn[iSamp]
            lmat[samp_coor[0],samp_coor[1],sampIdLab[iSamp]] += 1
        self.sampLab = sampLab
        self.sampIdLab = sampIdLab
        self.labToIdLab = labToIdLab
        self.idLabToLab = n.asarray(lab,dtype='O')

    def makeClassifierMatrix(self,classLabels=None):
        pass
        #if classLabels is None:
        #    classLabels = self.
