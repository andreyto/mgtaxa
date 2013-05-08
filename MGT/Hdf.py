### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##



import tables as pt

import os.path

def array_chunked_copy(src,dst,chunk):
    """Assign elements of src to elements of dst in chunks.
    PyTables cannot copy from Array or other Node directly like
    dst[:] = src. It requires intermediate conversion into 
    numpy array. If datasets are large, dst[:] = src[:]
    will use a lot of RAM. The method implemented here will be both fast
    and memory conserving of chunk parameter is selected properly.
    Alternative approach of using node.copy() is not appropriate when
    implicit element type conversion is expected to happen during 
    assignment.
    @param dst destination node
    @param source node
    @param chunk chunk size in records
    """
    l = len(src)
    for i in xrange(0,l,chunk):
        j = min(i+chunk,l)
        dst[i:j] = src[i:j]

def hdfSplitPath(path):
    return os.path.split(path)

def hdfDtype(hdfObj):
    """Return a numpy dtype corresponding to the PyTables datatype instance.
    Currently only tested on Table instances.
    @todo make sure it works for other datatypes."""
    return hdfObj.description._v_dtype

# Some PyTables gotchas:
# 1. Answer to the problem below is that after row.append(),
# the row object points to something past the end of table,
# and its content is undefined.
# If we have a Table.Row object and
#row['f'] = value
#then
#print row
#will print a tuple with old values, those before the asignment
#until we call table.flush() (row.append() alone does not help)
#row[:] print tuples of zeros even after table.flush() -
#slicing shows correct values only on pre-existing records.
#row.fetch_all_fields() behaves like 'print row'.
#row.nrow behaves equally.
#row['f'] always return the latest value
# 2. table.nrows is not updated until flush()
# 3. row.nrow is 0 before any row was appended (on empty table)
#and it will be 0 again after the first append()
# 4. Alternative to using table.row.append():
# a = numpy.zeros(1,dtype=table.description._v_dtype)
# table.append(a)
#

