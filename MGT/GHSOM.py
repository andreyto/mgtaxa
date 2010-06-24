### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Interface to GHSOM Hierarchical self-organizing maps program.
A Rauber, D Merkl, M Dittenbach - Neural Networks, IEEE Transactions on, 2002
The growing hierarchical self-organizing map: exploratory analysis of high-dimensional data.
http://www.ifs.tuwien.ac.at/~andi/ghsom/
"""

from MGT.Common import *
from MGT.SOM import *
import pdb

class GHSOM:
    def __init__(self,name):
        self.name = name
        self.ivFile = name+'.ghsom.iv'
        self.tvFile = name+'.ghsom.tv'
        self.mapIdStr = "_1_1_0_0"
        self.setModelDir(".")

    def setModelDir(self,dirName):
        makedir(dirName)
        self.modDir = dirName
    
    def getName(self,kind="unit"):
        return os.path.join(self.modDir,self.name+self.mapIdStr+"."+kind)

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

    def loadUnit(self):
        inp = open(self.getName(kind="unit"),'r')
        for line in inp:
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            parts = line.split(None,1)
            if len(parts) == 1:
                key = parts[0]
                dat = ""
            else:
                key, dat = parts
            assert key.startswith("$")
            key = key[1:]
            if key == 'TYPE':
                assert dat == "rect","Only rectangular grid currently supported"
            elif key == 'XDIM':
                xdim = int(dat)
                ydim = None
            elif key == 'YDIM':
                ydim = int(dat)
                grid = n.zeros((xdim,ydim),dtype='O')
                #for x in grid.flat:
                #    x = []
            elif key == 'POS_X':
                pos_x = int(dat)
                pos_y = None
            elif key == 'POS_Y':
                pos_y = int(dat)
            elif key == 'NR_VEC_MAPPED':
                nr_vec_mapped = int(dat)
            elif key == 'MAPPED_VECS':
                if nr_vec_mapped == 0:
                    mapped_vecs = n.zeros(nr_vec_mapped,dtype='O')
                else:
                    mapped_vecs = n.asarray(dat.split(),dtype='O')
                assert nr_vec_mapped == len(mapped_vecs)
                grid[pos_x,pos_y] = mapped_vecs
            elif key == 'NR_SOMS_MAPPED':
                assert int(dat) == 0, "Hierarchical SOMs are not yet supported"
                pos_x = None
                pos_y = None
                mapped_vecs = None
                nr_vec_mapped = None
        inp.close()
        return grid

    def loadWeights(self):
        inp = open(self.getName(kind="wgt"),'r')
        for line in inp:
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            parts = line.split(None,1)
            if len(parts) == 1:
                key = parts[0]
                dat = ""
            else:
                key, dat = parts
            assert key.startswith("$")
            key = key[1:]
            if key == 'TYPE':
                assert dat == "rect","Only rectangular grid currently supported"
            elif key == 'XDIM':
                xdim = int(dat)
                ydim = None
            elif key == 'YDIM':
                ydim = int(dat)
                grid = None
            elif key == 'VEC_DIM':
                vec_dim = int(dat)
                grid = n.zeros((xdim,ydim,vec_dim),dtype='f8')
                grid_filled = n.zeros((xdim,ydim),dtype='b')
                break
        for line in inp:
            if line == '\n' or line.startswith('#'):
                continue
            parts = line.rsplit(None,1)
            vec = n.fromstring(parts[0],sep=' ',dtype='f8')
            assert len(vec) == vec_dim
            # parse the row,col from "1_1_0/0_<row_number>/<col_number>"
            irow,icol = parts[1].split('/')[1:3]
            irow = int(irow.split('_')[1])
            icol = int(icol)
            grid[irow,icol] = vec
            assert grid_filled[irow,icol] == False
            grid_filled[irow,icol] = True
        assert grid_filled.all()
        inp.close()
        return grid

    def loadModel(self,components=("unit","weights")):
        mod = SOMModel()
        if "unit" in components:
            unit = self.loadUnit()
            mod.unit = unit
        if "weights" in components:
            weights = self.loadWeights()
            mod.weights = weights
        return mod


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
