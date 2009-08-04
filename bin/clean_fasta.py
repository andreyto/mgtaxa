### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.FastaIO import *


def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="append", type="string",dest="inSeq"),
        make_option("-l", "--line-len",
        action="store", type="int",dest="lineLen",default=30000),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

assert len(opt.inSeq) > 0, "Need at least one in-seq option"

for inSeq in opt.inSeq:
    inpSeq = FastaReader(inSeq)
    outSeq = open(inSeq+'.cl','w')
    iRec = 0
    for rec in inpSeq.records():
        hdr = rec.header()
        outSeq.write(hdr)
        for s in rec.seqChunks(opt.lineLen):
            outSeq.write(s)
            outSeq.write('\n')
        if iRec % 1000 == 0:
            print "Processed %i records from file %s" % (iRec,inSeq)
        iRec += 1

    outSeq.close()
    inpSeq.close()

