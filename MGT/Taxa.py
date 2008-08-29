"""Combines classes for processing taxonomy trees"""

from MGT.TaxaTree import *
from MGT.TaxaTreeDb import *
from MGT.TaxaIO import *
from MGT.TaxaIODb import *

def loadTaxaTree(ncbiDumpFile=options.taxaNodesFile,ncbiNamesDumpFile=options.taxaNamesFile):
    return TaxaTree(NodeStorageNcbiDump(ncbiDumpFile=ncbiDumpFile,
        ncbiNamesDumpFile=ncbiNamesDumpFile))

