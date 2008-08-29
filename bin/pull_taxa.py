"""Pull sequence for a given set of taxa directly from PANDA fasta file.
Created to test the correctness of MGTAXA own sequence/sample database.
The current implementation collects everything in memory."""

from MGT.Common import *
from itertools import izip
from MGT.Taxa import *

sampLen = 1000

pandaVirF = "/usr/local/db/panda/NuclSequences/ViralGroupNucl.fasta"

trSampF = "../taxaSamp.train.pkl"
tsSampF = "../taxaSamp.test.pkl"

labToIdF = "labelToId.pkl"

trSamp = loadObj(trSampF)
tsSamp = loadObj(tsSampF)

#giToTaxa = loadObj("/home/atovtchi/work/mgtdata/taxa.pkl.gz")

newTaxaDir="/home/atovtchi/work/mgtdata/taxonomy.new"

outFormat = "fasta"

#taxaTreeNew = loadTaxaTree(ncbiDumpFile=os.path.join(newTaxaDir,"nodes.dmp"),
#        ncbiNamesDumpFile=os.path.join(newTaxaDir,"names.dmp"))
#taxaTree = loadTaxaTree()

inpSeq = FastaReader(pandaVirF)

txSeq = {}

for samp in trSamp.values():
    for taxid in samp.keys():
        txSeq[taxid] = []

for samp in tsSamp.values():
    for taxid in samp.keys():
        txSeq[taxid] = []

iRec = 0
nTaxMis = 0
iRecSel = 0
nTaxMisSel = 0
for rec in inpSeq.records():
    hdr = rec.header()
    try:
        taxid = int(hdr.split("taxon:")[1].split(None,1)[0])
    except (ValueError,IndexError):
        taxid = 0
    gi = int(hdr.split("|",2)[2].split(None,1)[0])
    #if not giToTaxa[gi] == taxid:
    #    nTaxMis += 1
    if taxid in txSeq:
        txSeq[taxid].append(rec.sequence())
#        assert taxaTree.getNode(taxid).lineageRanksStr() == taxaTreeNew.getNode(taxid).lineageRanksStr()
#        if not giToTaxa[gi] == taxid:
#            nTaxMis += 1
        iRecSel += 1
    iRec += 1

print nTaxMis, " / ", iRec
print nTaxMisSel, " / ", iRecSel
inpSeq.close()

#for id in txSeq.keys():
#    txSeq[id] = 'NNNNNNN'.join(txSeq[id])

idLabAll = set()
for samp in (trSamp,tsSamp):
    for idsClass in samp.iterkeys():
        idLabAll.add(idsClass)
idLabAll = sorted(idLabAll)
if idLabAll[0] != 0:
    idLabAll.insert(0,0)
labToId = numpy.asarray(idLabAll,dtype='i4')

dumpObj(labToId,labToIdF)

#invert index to get id to label mapping

idToLab = {}
for (lab,id) in izip(xrange(1,len(labToId)),labToId[1:]):
    idToLab[id] = lab

if outFormat == "fasta":
    labForm = ">%i\n"
elif outFormat == "feature":
    labForm = "%i "
else:
    raise ValueError, "Unknown output format: " + outFormat

for outF,samp,outIdSampF in (("../train.p.svm",trSamp,"train.idsamp.pkl"),
        ("../test.p.svm",tsSamp,"test.idsamp.pkl")):
    out = open(outF,'w')
    idsamp = []
    for idClass in samp:
        lab = idToLab[idClass]
        for idSamp in samp[idClass]:
            for seq in txSeq[idSamp]:
                for sSamp in range(0,len(seq),sampLen):
                    eSamp = sSamp + sampLen
                    if eSamp <= len(seq):
                        out.write(labForm % lab)
                        out.write(seq[sSamp:eSamp])
                        out.write('\n')
                        idsamp.append(idSamp)
    out.close()
    idsamp = numpy.asarray(idsamp,dtype='i4')
    dumpObj(idsamp,outIdSampF)

