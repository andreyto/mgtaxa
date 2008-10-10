from MGT.GHSOM import *
from MGT.Svm import *
import pdb
import pylab as pl

som = GHSOM("all")
som.setModelDir(os.path.join(options.testDataDir,"som/wd_2_0/grow"))
mod = som.loadModel(components=("unit","weights"))
mod.makeUMatrix()
data = loadSparseSeqsAsDense(inpFile=os.path.join(options.testDataDir,"som/wd_2_0/all.svm"))
mod.setSamples(data)
print "Mapping samples..."
mod.mapSamples(pmRadius=mod.paretoRadiusLenNormSamp)
#mod.mapSamples(pmRadius=0.5)
mod.makeUStarMatrix()
pl.imshow(mod.pmat.T,aspect="equal",interpolation="nearest")
pl.savefig("tmp.png")
#mod.paretoRadius()
