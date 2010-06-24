### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Split a multi-FASTA file into chunks constrained by a total size"""

__all__ = [ "splitFasta" ]

from MGT.FastaIO import *
from MGT.Common import*
#from MGT.Shogun import *
#from MGT.Shogun.Util import *
#from shogun.Features import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-fasta",
        action="store", type="string",dest="inpFasta"),
        make_option("-o", "--out-name",
        action="store", type="string",dest="outName"),
        make_option("-n", "--chunk-size",
        action="store", type="string",dest="chunkSize",default="500M"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


def splitFasta(inpFasta,outName,chunkSize):
    inpSeq = FastaReader(inpFasta)
    iChunk = 1
    out = openCompressed("%s-%05i" % (outName,iChunk),"w")
    print "Writing chunk %i of target size %i" % (iChunk,chunkSize)
    for rec in inpSeq.records():
        hdr = rec.header()
        out.write(hdr)
        lenSeq = 0
        for chunk in rec.seqChunks(32*1024):
            out.write(chunk)
            out.write("\n")
            lenSeq += len(chunk)
        # we approximate the next seq length by the last one
        if out.tell() + lenSeq >= chunkSize:
            out.close()
            iChunk += 1
            out = openCompressed("%s-%05i" % (outName,iChunk),"w")
            print "Writing chunk %i of target size %i" % (iChunk,chunkSize)
    inpSeq.close()


if __name__ == "__main__":

    opt,args = getProgOptions()
    assert opt.inpFasta is not None and opt.outName is not None

    chunkSize = opt.chunkSize.upper()
    if chunkSize.endswith("M"):
        chunkSize = int(chunkSize[:-1])*(10**6)
    elif chunkSize.endswith("G"):
        chunkSize = int(chunkSize[:-1])*(10**9)
    elif chunkSize.endswith("K"):
        chunkSize = int(chunkSize[:-1])*(10**3)
    elif not chunkSize[-1].isdigit():
        raise ValueError("--chunk-size %s has illegal suffix" % chunkSize)
    else:
        chunkSize = int(chunkSize)
    splitFasta(inpFasta=opt.inpFasta,outName=opt.outName,chunkSize=chunkSize)

