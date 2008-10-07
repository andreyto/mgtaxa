from MGT.GHSOM import *
import pdb

som = GHSOM("all")
som.setModelDir(os.path.join(options.testDataDir,"som/wd_2_0/grow"))
mod = som.loadModel(components=("unit","weights"))
mod.makeUMatrix()
pdb.set_trace()

