### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *

a = numpy.arange(10,dtype=int)
b = numpy.ones(10,dtype=int)

a[-1] = 15
a[-2] = 16

b[2] = 10
b[5] = 15

c = fromWhereItems(whItems = {'ind':a,'val':b},defVal=0)

for i,v in zip(a,b):
    assert c[i] == v

pdb.set_trace()
