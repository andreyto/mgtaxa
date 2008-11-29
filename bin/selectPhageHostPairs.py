"""Select phage/host RefSeq record pairs"""
from MGT.Taxa import *
from MGT.FastaIO import *
from MGT.Phage import *
import pdb


def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="append", type="string",dest="inSeq"),
        make_option("-p", "--hosts",
        action="store", type="string",dest="inHosts"),
        make_option("-o", "--out-picks",
        action="store", type="string",dest="outPicks"),
        make_option("-s", "--out-seq-ids",
        action="store", type="string",dest="outSeqIds"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()


giToTaxa = loadGiTaxBinNew()
taxaTree = loadTaxaTreeNew()


#taxaTree.getRootNode().setIsUnderUnclass()

mapFastaRecordsToTaxaTree(inSeqs=opt.inSeq,taxaTree=taxaTree,giToTaxa=giToTaxa)
loadHosts(inpFile=opt.inHosts,taxaTree=taxaTree)

seqPicker = PhageHostSeqPicker(taxaTree=taxaTree)
seqPicker.pickSeqHosts()
seqPicker.groupSeqHosts()
seqPicker.printGroupSeqHosts()

seqPicker.pickPairs(maxMicSpe=1,maxMicSeq=1,maxVir=1,maxVirSeq=1)
seqPicker.checkPairs(giToTaxa=giToTaxa)

seqPicker.save(opt.outPicks)
seqPicker.load(opt.outPicks)
seqPicker.checkPairs(giToTaxa=giToTaxa)
seqPicker.saveSeqIds(opt.outSeqIds)


