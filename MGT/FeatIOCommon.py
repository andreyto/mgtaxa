### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.FeatCommon import *

def featIdFileNameDef(featFile):
    return featFile+'.id'


def svmSaveId(ids,out):
    if not hasattr(out,'write'):
        out = openCompressed(out,'w')
        ownOut = True
    else:
        ownOut = False
    s = ''.join([ "%s\n" % x for x in ids ])
    out.write(s)
    if ownOut:
        out.close()

def svmLoadId(inp,dtype=idDtype):
    if not hasattr(inp,'read'):
        inp = openCompressed(inp,'r')
        ownInp = True
    else:
        ownInp = True
    x = n.asarray([ s.strip() for s in inp ],dtype=dtype)
    if ownInp:
        inp.close()
    return x
    
