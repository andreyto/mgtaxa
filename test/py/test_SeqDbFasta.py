from MGT.SeqDbFasta import *

seqDbPath = pjoin(options.testDataDir,"seqdb-fasta")

seqDb = SeqDbFasta(seqDbPath)
ids = seqDb.getIdList()
for id in ids:
    seqDb.writeFastaBothStrands([id],"%s.tmp.fna.gz" % (id,))


