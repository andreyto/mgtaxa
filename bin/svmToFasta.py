### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Svm import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="store", type="string",dest="inSeq"),
        make_option("-o", "--out-seq",
        action="store", type="string",dest="outSeq"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

data = loadSeqs(opt.inSeq)
out = openCompressed(opt.outSeq,'w')

for rec in data:
    out.write(">%s\n%s\n" % (rec['id'],rec['feature']))
out.close()

