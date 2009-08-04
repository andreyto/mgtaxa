### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import tables
import numpy

fileh = tables.openFile('earray1.h5', mode='w')
a = tables.StringAtom(itemsize=1)
# Use ``a`` as the object type for the enlargeable array.
array_c = fileh.createEArray(fileh.root,
        'array_c',
        a,
        (0,),
        "Chars",
        expectedrows=1000,
        filters=tables.Filters(complevel=1, complib='lzo',shuffle=False)
        )
#size of internal memory buffer of this leaf
array_c.nrowsinbuf = 1000
#array_c.append(numpy.array(['a'*2, 'b'*4], dtype='S8'))
n_array = numpy.fromstring('a'*10,dtype='S1')
print n_array
array_c.append(n_array)
#This will not work - it things that I am supplying a single atom of incorrect size
#array_c.append('b'*20)
array_c.append(['b']*20)

# Read the string ``EArray`` we have created on disk.
#for s in array_c:
#    print 'array_c[%s] => %r' % (array_c.nrow, s)
print array_c[2:7]
# Close the file.
fileh.close()
