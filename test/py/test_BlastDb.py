### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.BlastDb import *

from MGT.Util import strToFile
from MGT.Config import options

gis = """
62184494
58430915
"""
giFile = "tmp_gi.list"
dbAlias = "tmp_test"
strToFile(gis,giFile)

blastDb = BlastDb()
#blastDb.makeDbAlias(dbNameAlias=dbAlias,giFile=giFile)

#inp = blastDb.fastaReader(dbName=dbAlias,giFile=giFile)
#for rec in inp.records():
    #print rec.header(),
#inp.close()

blastDb.makeDbAlias(dbNameAlias=dbAlias,giFile=None)

inp = blastDb.fastaReader(dbName=dbAlias,giFile=None)
for rec in inp.records():
    hdr = rec.header()
    #seqlen = rec.seqLen()
    seqlen = 0
    for chunk in rec.seqArrays(chunkSize=100000):
        seqlen += len(chunk)
    print hdr.rstrip("\n") + " len:%s" % (seqlen,)
inp.close()
