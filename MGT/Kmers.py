from MGT.Common import *

try:
    from MGTX.kmersx import *
except ImportError:
    class SvmSparseFeatureWriterTxt:
        pass

