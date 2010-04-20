from MGT.SeqDbFasta import *

seqDbPath = pjoin(options.testDataDir,"seqdb-fasta")

seqDb = SeqDbFasta(seqDbPath)
ids = seqDb.getIdList()
for id in ids:
    seqDb.writeFastaBothStrains([id],"%s.tmp.fna.gz" % (id,))


