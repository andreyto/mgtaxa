from MGT.FastaIO import *
from MGT.SeqRandomApp import *

inpFile = '/usr/local/projects/GOSII/shannon/Indian_Ocean_Viral/asm_combined_454_large/454LargeContigs.fna'
wr = FastaWriter("asm_combined_454_large.5K.rnd.fna",lineLen=80)
#this will only shuffle real sequences, maintaining 1-mer frequencies (primarily GC content, in other words)
rndgen = RandomSeqGenShuffle()
for rec in FastaReader(open(inpFile,'r')).records():
    hdr = rec.header()
    seq = rec.sequence()
    if len(seq) >= 5000:
        wr.record(hdr,rndgen(seq).tostring())
wr.close()

