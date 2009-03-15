"""Create binary representation of NCBI taxonomy dump."""

from MGT.Taxa import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-dump",
        action="append", type="string",dest="inDump"),
        make_option("-o", "--out-bin-file",
        action="store", type="string",dest="outFile"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()


makeGiTaxBin(opt.inDump,opt.outFile)

