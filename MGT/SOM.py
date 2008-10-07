"""Operations on Self Organizing Maps.
The Maps are trained by specific implementation defined elsewhere.
Here we implement some common post-processing tasks such as
U-matrix computation, segmentation and clustering on the maps."""

from MGT.Common import *

class SOMModel(Struct):
    
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
        u = n.zeros(w.shape[:2],dtype='f8')
        c = n.zeros(w.shape[:2],dtype='i4')
        #dx = n.zeros((w.shape[0]-1,w.shape[1]  ),dtype='f8')
        #dy = n.zeros((w.shape[0]  ,w.shape[1]-1),dtype='f8')
        #dd = n.zeros((w.shape[0]-1,w.shape[1]-1),dtype='f8')
        dx = n.sqrt(((w[1:, :,:] - w[:-1,:  ,:])**2).sum(axis=2))
        dy = n.sqrt(((w[: ,1:,:] - w[:  ,:-1,:])**2).sum(axis=2))
        dd = n.sqrt(((w[1:,1:,:] - w[:-1,:-1,:])**2).sum(axis=2))
        u[:-1,:] += dx
        c[:-1,:] += 1
        u[1:,:] += dx
        c[1:,:] += 1
        u[:,:-1] += dy
        c[:,:-1] += 1
        u[:,1:] += dy
        c[:,1:] += 1
        u[:-1,:-1] += dd
        c[:-1,:-1] += 1
        u[1:,1:] += dd
        c[1:,1:] += 1
        u[c>0] /= c[c>0]
        self.umat = u

