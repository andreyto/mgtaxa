from MGT.Common import *
import pdb

"""Interface to GHSOM Hierarchical self-organizing maps program.
A Rauber, D Merkl, M Dittenbach - Neural Networks, IEEE Transactions on, 2002
The growing hierarchical self-organizing map: exploratory analysis of high-dimensional data.
http://www.ifs.tuwien.ac.at/~andi/ghsom/
"""

class GHSOM:
    def __init__(self,name):
        self.name = name
        self.ivFile = name+'.ghsom.iv'
        self.tvFile = name+'.ghsom.tv'

    def writeInput(self,data):
        m = data['feature']
        l = data['label']
        nVec = m.shape[0]
        nFeat = m.shape[1]
        ivout = open(self.ivFile,'w')
        ivout.write(ivHdr % dict(nVec=nVec,nFeat=nFeat))
        frmt = "%f "*nFeat + "%s\n"
        for iVec in xrange(nVec):
            ivout.write(frmt % (tuple(m[iVec])+(l[iVec],)))
        ivout.close()
        # Each row is
        #0 <name of feature 1> <df> <tf> <min_tf> <max_tf> <mean_tf>
        #1 .....
        # We just fill with default values
        tvm = numpy.ones((nFeat,7),dtype='i4')
        tvm[:,0] = numpy.arange(nFeat,dtype='i4')
        tvm[:,1] = tvm[:,0]
        s = "\n".join([ " ".join([ "%s" % x for x in row ]) for row in tvm ])
        tvout = open(self.tvFile,'w')
        tvout.write(tvHdr % dict(nVec=nVec,nFeat=nFeat))
        tvout.write(s)
        tvout.close()

    def writeInputFromShogunSparse(self,feat,lab):
        f = feat.get_full_feature_matrix().transpose()
        self.writeInput(data=dict(feature=f,label=lab))


ivHdr=\
"""$TYPE inputvec
$XDIM %(nVec)i
$YDIM 1
$VECDIM %(nFeat)i
"""

tvHdr=\
"""$TYPE template
$XDIM 7
$YDIM  %(nVec)i
$VECDIM %(nFeat)i
"""
