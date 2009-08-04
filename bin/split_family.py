### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Pull sequence for a given set of taxa directly from PANDA fasta file.
Created to test the correctness of MGTAXA own sequence/sample database.
The current implementation collects everything in memory."""

from MGT.Common import *
from itertools import izip
from MGT.Taxa import *

pandaVirF = "/usr/local/db/panda/NuclSequences/ViralGroupNucl.fasta"

trSampF = "../taxaSamp.train.pkl"
tsSampF = "../taxaSamp.test.pkl"

labToIdF = "labelToId.pkl"

txSeqF = "txSeq.pkl"

trSamp = loadObj(trSampF)
tsSamp = loadObj(tsSampF)

#giToTaxa = loadObj("/home/atovtchi/work/mgtdata/taxa.pkl.gz")

newTaxaDir="/home/atovtchi/work/mgtdata/taxonomy.new"

#taxaTreeNew = loadTaxaTree(ncbiDumpFile=os.path.join(newTaxaDir,"nodes.dmp"),
#        ncbiNamesDumpFile=os.path.join(newTaxaDir,"names.dmp"))
taxaTree = loadTaxaTree()

if os.path.isfile(txSeqF):

    txSeq = loadObj(txSeqF)

else:

    inpSeq = FastaReader(pandaVirF)

    txSeq = {}

    for samp in trSamp.values():
        for taxid in samp.keys():
            txSeq[taxid] = []

    for samp in tsSamp.values():
        for taxid in samp.keys():
            txSeq[taxid] = []

    iRec = 0
    iRecSel = 0
    for rec in inpSeq.records():
        hdr = rec.header()
        try:
            taxid = int(hdr.split("taxon:")[1].split(None,1)[0])
        except (ValueError,IndexError):
            taxid = 0
        #gi = int(hdr.split("|",2)[2].split(None,1)[0])
        if taxid in txSeq:
            txSeq[taxid].append(rec.sequence())
            iRecSel += 1
        iRec += 1

    print iRec
    print iRecSel
    inpSeq.close()

    for id in txSeq.keys():
        txSeq[id] = 'NNNNNNN'.join(txSeq[id])

    dumpObj(txSeq,txSeqF)

# On fam1:fam2 we got 
#
#cm      =
#        m       =       [[   0    0    0]
#                         [   0 1312   11]
#                         [   0  359  148]]
#        mb      =       [[   0    0    0]
#                         [   0 1312   11]
#                         [   0  718  296]]
#
#sen     =       [ 0.99168556  0.29191321]
#senLab  =       [1 2]
#senMean =       0.641799389052
#spe     =       [ 0.64630542  0.96416938]
#speLab  =       [1 2]
#speMean =       0.805237399913
#
#[lab:id:spe...]:  1:10780:0.65  2:10811:0.96
#
#[lab:id:sen...]:  1:10780:0.99  2:10811:0.29
#
# 10780 attracts both own and the other's 
# samples.
# So, we split 10780 into genus labels and
# keep 10811 as it is. We save the sequence taxid
# for each sample as well. That will allow us to
# see how the misassigned sample are distributed.

spId = 10780
kpId = 10811

spSamp = (trSamp[spId],tsSamp[spId])

newSamp = ({kpId:trSamp[kpId]},{kpId:tsSamp[kpId]})

newRank = "genus"

for n,sp in zip(newSamp,spSamp):
    for id in sp.keys():
        newIdNode = taxaTree.getNode(id).findRankInLineage(rank=newRank)
        if newIdNode is not None:
            newId = newIdNode.id
            if not n.has_key(newId):
                n[newId] = {}
            n[newId][id] = sp[id]

idLabAll = set()
for samp in newSamp:
    for idsClass in samp.iterkeys():
        idLabAll.add(idsClass)
idLabAll= sorted(idLabAll)
if idLabAll[0] != 0:
    idLabAll.insert(0,0)
labToId = numpy.asarray(idLabAll,dtype='i4')

dumpObj(labToId,labToIdF)

#invert index to get id to label mapping

idToLab = {}
for (lab,id) in izip(xrange(1,len(labToId)),labToId[1:]):
    idToLab[id] = lab

for outF,samp,outIdSampF in (("../train.p.svm",newSamp[0],"train.idsamp.pkl"),
        ("../test.p.svm",newSamp[1],"test.idsamp.pkl")):
    out = open(outF,'w')
    idsamp = []
    for idClass in sorted(samp.keys()):
        lab = idToLab[idClass]
        for idSamp in sorted(samp[idClass].keys()):
            if len(txSeq[idSamp]) > 0:
                out.write(">%i\n" % lab)
                out.write(txSeq[idSamp])
                out.write('\n')
                idsamp.append(idSamp)
    out.close()
    idsamp = numpy.asarray(idsamp,dtype='i4')
    dumpObj(idsamp,outIdSampF)

