### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Shred each a record in a multi-FASTA file into fragments of a given size"""


from MGT.FastaIO import *
from MGT.Common import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-fasta",
        action="store", type="string",dest="inpFasta"),
        make_option("-o", "--out-fasta",
        action="store", type="string",dest="outFasta"),
        make_option(None, "--frag-size",
        action="store", type="int",dest="fragSize",default=400),
        make_option(None, "--frag-count-ratio",
        action="store", type="float",dest="fragCountRatio",default=1.),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


if __name__ == "__main__":

    opt,args = getProgOptions()
    assert opt.inpFasta is not None and opt.outFasta is not None

    shredFasta(inpFasta=opt.inpFasta,outFasta=opt.outFasta,
            fragSize=opt.fragSize,fragCountRatio=opt.fragCountRatio)

