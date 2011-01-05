from MGT.FastaIO import *

hdrSel=">scf1118659309935"
inpFile = "scaff-gos-bac/nonViral40kbPlusScaffolds.fa.gz"
wr = FastaWriter("scf1118659309935.fna",lineLen=80)
for rec in FastaReader(inpFile).records():
    hdr = rec.header().strip()
    seq = rec.sequence()
    if hdr==hdrSel:
        wr.record(hdr,seq)
wr.close()

