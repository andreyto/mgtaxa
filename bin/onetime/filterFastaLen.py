from MGT.FastaIO import *

inpFile = '/usr/local/projects/GOSII/shannon/Indian_Ocean_Viral/asm_combined_454_large/454LargeContigs.fna'
wr = FastaWriter("asm_combined_454_large.5K.fna",lineLen=80)
for rec in FastaReader(open(inpFile,'r')).records():
    hdr = rec.header()
    seq = rec.sequence()
    if len(seq) >= 5000:
        wr.record(hdr,seq)
wr.close()

