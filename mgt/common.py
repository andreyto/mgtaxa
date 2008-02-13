from AS_SCM_util import *

import numpy
import numpy.random as nrnd
from numpy import array
import os, sys, atexit, re, gzip
import itertools

import pdb

class PhyOptions:
    def __init__(self):
        self.taxaPickled = 'taxa.pkl'
        self.fastaHdrSqlLen = 40
        self.blastDataDir = 'blast'
        self.selGiFile = 'phyla_sel.gi'
        self.selDumpFile = 'phyla_sel.csv'
        self.selFastaFile = 'phyla_sel.fasta.gz'
        self.srcDbNameAlias = 'phyla'
        self.taxaDataDir = 'taxonomy'
        self.taxaCatFile = os.path.join(self.taxaDataDir,'taxcat','categories.dmp')
        self.taxaGiFile = os.path.join(self.taxaDataDir,'gi_taxid_nucl.dmp') 
        self.taxaDumpDir = os.path.join(self.taxaDataDir,'taxdump')
        self.taxaNodesFile = os.path.join(self.taxaDumpDir,'nodes.dmp')
        self.taxaDivisionFile = os.path.join(self.taxaDumpDir,'division.dmp')
        self.kmerTxtDir = os.environ['PHYLA_KMERS']
        self.kmerTestFile = os.path.join(self.kmerTxtDir,'6mers_1K.gz')


def masksToInd(maskVal,maskInd):
    """Convert value[0:someN] and index[0:someN] into value[index]. 
    Assumes N->1 relation between index and value. In particular, that means
    value[i] == 0 && value[j] != 0 && index[i] == index[j] is not allowed.
    No checking is made for the above prerequsite - the caller is responsible."""
    val = numpy.zeros(numpy.max(maskInd)+1,dtype=maskVal.dtype)
    val[maskInd] = maskVal
    return val
    
def whereItems(arr,condition):
    wh = numpy.where(condition)
    return numpy.rec.fromarrays((wh[0],arr[wh]),names='ind,val')

def taxidFromPhyFastaHeader(hdr):
    """Extract taxid from a FASTA header formatted for PHY database"""
    return int(hdr.split("taxid:",1)[1].split(" ",1)[0])

def joinFastaByTaxid(inpFastaFile,outFastaFile):
    reader = fastaReaderGzip(inpFastaFile)
    writer = openGzip(outFastaFile,'w')
    taxidLast = -1
    for rec in reader.records():
        taxid = taxidFromPhyFastaHeader(rec.header())
        if taxid != taxidLast:
            writer.write(">%s\n" % (taxid,))
            taxidLast = taxid
        for line in rec.seqLines():
            writer.write(line)
    reader.close()
    writer.close()
