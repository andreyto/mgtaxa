
import tables as pt

import os.path

def hdfSplitPath(path):
    return os.path.split(path)

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

