### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Create binary representation of NCBI taxonomy dump."""

from MGT.Taxa import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-o", "--out-bin-file",
        action="store", type="string",dest="outFile"),
    ]
    parser = OptionParser(usage = \
            "usage: %prog [options] gi-to-taxa.dmp [,gi-to-taxa.dmp...]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

makeGiTaxBin(args,opt.outFile)

