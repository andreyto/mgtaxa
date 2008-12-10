"""Select phage/host RefSeq record pairs"""
from MGT.SeqIO import *
import pdb


def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="append", type="string",dest="inSeq"),
        make_option("-o", "--out-seq",
        action="store", type="string",dest="outSeq"),
        make_option("-s", "--in-seq-ids",
        action="store", type="string",dest="inSeqIds"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()


#taxaTree = loadTaxaTreeNew()

seqIds = loadObj(opt.inSeqIds)

pullNCBISeq(inSeqs=opt.inSeq,seqIds=seqIds,outSeq=opt.outSeq)

