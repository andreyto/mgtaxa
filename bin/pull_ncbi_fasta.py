### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Read a file with sequence features and create a file with sparse real features based on word distance histogram."""
from MGT.Taxa import *
from MGT.Svm import *
from MGT.FastaIO import *
from MGT.Shogun.Util import *
from shogun.Features import *
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="append", type="string",dest="inSeq"),
        make_option("-o", "--out-name",
        action="store", type="string",dest="outName"),
        make_option("-s", "--num-splits",
        action="store", type="int",dest="nSplits",default=3),
        make_option("-m", "--min-samp-count",
        action="store", type="int",dest="minSampCount",default=100),
        make_option("-t", "--max-samp-seq",
        action="store", type="int",dest="maxSampCountPerSeq"),
        make_option("-l", "--samp-len",
        action="store", type="int",dest="sampLen",default=1000),
        make_option("-f", "--samp-offset",
        action="store", type="int",dest="sampOffset",default=0),
        make_option("-d", "--make-other",
        action="store_true", dest="makeOther",default=False),
        make_option("-a", "--alphabet",
        action="store", type="choice",choices=("dna","protein"),
        dest="alphabet",default="dna"),
        make_option("-e", "--degen-len",
        action="store", type="int",dest="degenLen",default=1),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

if opt.alphabet == "dna" and opt.degenLen >= 0:
    nonDegenSymb = 'ACGT'
    symCompr = SymbolRunsCompressor('N',opt.degenLen)
else:
    symCompr = lambda s: s

nSplits = opt.nSplits

doWriteOtherGroup = opt.makeOther

mgtDbDir = "/home/atovtchi/work/mgtdata"

taxaDir=os.path.join(mgtDbDir,"taxonomy.new")

giToTaxa = loadGiTaxBin(os.path.join(taxaDir,"gi_taxid.pkl.gz"))

taxaTree = loadTaxaTree(ncbiDumpFile=os.path.join(taxaDir,"nodes.dmp"),
        ncbiNamesDumpFile=os.path.join(taxaDir,"names.dmp"))
#taxaTree = loadTaxaTree()

taxaTree.getRootNode().setIsUnderUnclass()
viralNode = taxaTree.getNode(viralRootTaxid)
viroidsNode = taxaTree.getNode(12884)
dsDnaNode = taxaTree.getNode(35237)
cellNode = taxaTree.getNode(131567) # cellular organisms - bact, arch & euk
phageTailedNode = taxaTree.getNode(phageTailedTaxid)
#the order is important here because phageTailedNode is a subnode of viralNode
#whichSuperNode should find phages first
famNodes = (phageTailedNode,viralNode,cellNode)

topNode = phageTailedNode

#viralTaxidLev2 = (\
#        35237, #dsDNA
#        35325, #dsRNA
#        35268, #retroid
#        29258, #ssDNA
#        439488, #ssRNA
#        )

viralNodesLev2 = [ taxaTree.getNode(id) for id in viralTaxidLev2 ]

def whichSupernode(node,others):
    for o in others:
        if node.isSubnode(o):
            return o
    return None


otherGroup = Struct()
otherGroup.splitSeq = numpy.zeros(nSplits,'O')
for i in range(nSplits):
    otherGroup.splitSeq[i] = []

iRec = 0
taxMis = Struct(bounds=0,zeroG=0,trN=0,trV=0)
for inSeq in opt.inSeq:
    inpSeq = FastaReader(inSeq)
    for rec in inpSeq.records():
        hdr = rec.header()
        gi = rec.getNCBI_GI()
        if len(giToTaxa) <= gi:
            taxMis.bounds += 1
            print "giToTaxa bounds: "+hdr
        else:
            taxid = giToTaxa[gi]
            if taxid == 0:
                taxMis.zeroG += 1
                print "zero giToTaxa: "+hdr
            else:
                try:
                    node = taxaTree.getNode(taxid)
                except KeyError:
                    taxMis.trN += 1
                    print "no node %s %s" % (taxid,hdr)
                else:
                    if False and (not node.isSubnode(viralNode) and not node.isSubnode(viroidsNode)):
                        taxMis.trV += 1
                        print "not under viral node %s %s %s" % (taxid,node,hdr)
                    else:
                        if node.isSubnode(topNode):
                            if not node.isUnderUnclass:
                                genNode = node.findRankInLineage("species")
                                if genNode is not None:
                                    famNode = genNode.findRankInLineage("genus")
                                    #famNode = whichSupernode(genNode,famNodes) #viralNodesLev2
                                    if famNode is not None:
                                        if not hasattr(famNode,'splitCounts'):
                                            famNode.splitCounts = numpy.zeros(nSplits,'i4')
                                            famNode.splitSeq = numpy.zeros(nSplits,'O')
                                            for i in range(nSplits):
                                                famNode.splitSeq[i] = []
                                        try:
                                            splitId = genNode.splitId
                                        except AttributeError:
                                            scnt = famNode.splitCounts
                                            # among equal minimal counts, add to the randomly picked
                                            splitIdMin = scnt.argmin()
                                            splitIdsMin = numpy.where(scnt==scnt[splitIdMin])[0]
                                            splitId = nrnd.permutation(splitIdsMin)[0]
                                            genNode.splitId = splitId
                                        famNode.splitCounts[splitId]+=1
                                        seq = symCompr(rec.sequence())
                                        stride = opt.sampLen + opt.sampOffset
                                        sampStarts = nrnd.permutation(range(0,len(seq)-opt.sampLen,stride))
                                        if opt.maxSampCountPerSeq < len(sampStarts):
                                            sampStarts = sampStarts[:opt.maxSampCountPerSeq]
                                        for sampStart in sampStarts:
                                            samp = Struct(seq=seq[sampStart:sampStart+opt.sampLen],
                                                    id=node.id)
                                            if checkSaneAlphaHist(samp.seq,nonDegenSymb,minNonDegenRatio=0.9):
                                            #if len(samp) >= opt.sampLen*0.9:
                                                famNode.splitSeq[splitId].append(samp)
                                            else:
                                                print "Too many degenerate symbols in sample, skipping: ",samp
        iRec += 1
        if iRec % 1000 == 0:
            print "Scanned %s records" % iRec
    inpSeq.close()

print "Scanned %d records. Taxonomy mismatches:\n%s\n" % (iRec,taxMis)

print "Writing splits..."

svmWriters = [SvmStringFeatureWriterTxt(opt.outName+'.s.%d' % iSplit) 
        for iSplit in xrange(nSplits)]
if doWriteOtherGroup:
    labelOther = 1
else:
    labelOther = 0
label = labelOther + 1

labToId = numpy.zeros(10000,'i4')
for node in topNode.iterDepthTop():
    if hasattr(node,'splitCounts'):
        isOther = True
        if (node.splitCounts > 0).all():
            sampCounts = numpy.asarray([ len(l) for l in node.splitSeq ])
            if (sampCounts >= opt.minSampCount).all():
                try:
                    labToId[label] = node.id
                except:
                    labToId.resize(label*2)
                    labToId[label] = node.id
                print "Writing label %i for node %s" % (label,node.lineageStr())
                for iSplit in xrange(nSplits):
                    print "%s samples in split %s for label %s" % \
                            (len(node.splitSeq[iSplit]),iSplit,label)
                    for samp in node.splitSeq[iSplit]:
                       svmWriters[iSplit].write(label,samp.seq,samp.id)
                label += 1
                isOther = False
        if isOther:
            iSplitOther = nrnd.randint(nSplits)
            for iSplit in xrange(nSplits):
                otherGroup.splitSeq[iSplitOther].extend(node.splitSeq[iSplit])

if doWriteOtherGroup:
    print "Writing label %i for node Other" % (labelOther,)
    for iSplit in xrange(nSplits):
        print "%s samples in split %s for label %s" % \
                (len(otherGroup.splitSeq[iSplit]),iSplit,labelOther)
        for samp in otherGroup.splitSeq[iSplit]:
           svmWriters[iSplit].write(labelOther,samp.seq,samp.id)

for svmWriter in svmWriters:
    svmWriter.close()

labToId.resize(label+1)
dumpObj(labToId,opt.outName+'.labid')

print "Labels-TaxIds: ", labToId

