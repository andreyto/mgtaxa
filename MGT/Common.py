from MGT.Util import *
from MGT.Config import Options, options

import numpy
import numpy.random as nrnd
from numpy import array
from random import sample as sampleWOR #Random sampling w/o replacement
import os, sys, atexit, re, gzip
import itertools

import pdb

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

def fromWhereItems(whItems,defVal=0):
    wh =  whItems['ind']
    a = numpy.empty(numpy.max(wh) + 1, dtype = whItems['val'].dtype)
    a[:] = defVal
    a[wh] = whItems['val']
    return a


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
