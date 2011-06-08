from MGT.FastaIO import *
import re
import sys

hdrSel=sys.argv[1]
inpFile = sys.stdin
wr = FastaWriter(sys.stdout,lineLen=80)
for rec in FastaReader(inpFile).records():
    hdr = rec.header().strip()
    if re.match(hdrSel,hdr[1:]):
        seq = rec.sequence()
        wr.record(hdr,seq)
#wr.close()

