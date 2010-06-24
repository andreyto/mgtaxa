### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from AS_SCM_util import *

from StringIO import StringIO


fastaInput = \
""">gi|152066|gb|M35122.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40987 seq_len:1040 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHAFIXA A.caulinodans nitrogen fixation protein (nifO and fixA) genes, complete cds and 5'end
GAGCTCGGCCTCTATGACATCGACGCCAGCGCGGTGAACGTCGCGCACGTGCCCGTCATTCCGGACGAGAACGAGGTGAG
GGCGGAGCCGGGGCGGAAAGCGCATTGTGGCATGCCAGACAGCCCTTTGATTTCATGCGCGTTTTCGGGCTGAAAGACAG
TTGGTACGACACTTGCTCATTCCTCCCCAAGAGCCCAACCGTTCCGGGAGCGAACGCAATGCACATCGTCGTCTGCATCA
AGCAGGTTCCTGACTCCGCGCAGATCCGCGTGCACCCCGTGACGAACACCATCATGCGTCAGGGTGTGCCCACGATCATC
>gi|152069|gb|M60872.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40988 seq_len:1110 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHANODDA A.caulinodans nodD gene, complete cds
TCGTGCAGAGATACCATGCGCTGTGCGCCTACAGCGAAGCAAGATGTGACGCTGGACAATCTTTCGTAGCTACGGCAAAA
TTCACCGTGGGACGACACTGCGGATGCTGGGTAGATGGAGAGTTAGTTGCGATTTAAGGGACTTGATCTGAATCTGCTTG
CAGCGCAGAATATTCGGGAGCTTCCATATAGTCGTTTTCCGATATGCAGAGATCAAGAGCTCTCAATCGC
>gi|443814|emb|X77126.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:262282 seq_len:982 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 AC16S6829 A.caulinodans (LMG 6829) 16S rRNA gene
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGCTCTTTCGCCGGTGAAGATAATGACGGT
GGTGCTGCATGGCTGTCGTCAG
"""

outExpected = \
"""
MODE = hdr *******************************************************************

>gi|152066|gb|M35122.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40987 seq_len:1040 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHAFIXA A.caulinodans nitrogen fixation protein (nifO and fixA) genes, complete cds and 5'end
>gi|152069|gb|M60872.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40988 seq_len:1110 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHANODDA A.caulinodans nodD gene, complete cds
>gi|443814|emb|X77126.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:262282 seq_len:982 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 AC16S6829 A.caulinodans (LMG 6829) 16S rRNA gene

MODE = seq *******************************************************************

>gi|152066|gb|M35122.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40987 seq_len:1040 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHAFIXA A.caulinodans nitrogen fixation protein (nifO and fixA) genes, complete cds and 5'end
GAGCTCGGCCTCTATGACATCGACGCCAGCGCGGTGAACGTCGCGCACGTGCCCGTCATTCCGGACGAGAACGAGGTGAG
GGCGGAGCCGGGGCGGAAAGCGCATTGTGGCATGCCAGACAGCCCTTTGATTTCATGCGCGTTTTCGGGCTGAAAGACAG
TTGGTACGACACTTGCTCATTCCTCCCCAAGAGCCCAACCGTTCCGGGAGCGAACGCAATGCACATCGTCGTCTGCATCA
AGCAGGTTCCTGACTCCGCGCAGATCCGCGTGCACCCCGTGACGAACACCATCATGCGTCAGGGTGTGCCCACGATCATC
>gi|152069|gb|M60872.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40988 seq_len:1110 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHANODDA A.caulinodans nodD gene, complete cds
TCGTGCAGAGATACCATGCGCTGTGCGCCTACAGCGAAGCAAGATGTGACGCTGGACAATCTTTCGTAGCTACGGCAAAA
TTCACCGTGGGACGACACTGCGGATGCTGGGTAGATGGAGAGTTAGTTGCGATTTAAGGGACTTGATCTGAATCTGCTTG
CAGCGCAGAATATTCGGGAGCTTCCATATAGTCGTTTTCCGATATGCAGAGATCAAGAGCTCTCAATCGC
>gi|443814|emb|X77126.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:262282 seq_len:982 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 AC16S6829 A.caulinodans (LMG 6829) 16S rRNA gene
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGCTCTTTCGCCGGTGAAGATAATGACGGT
GGTGCTGCATGGCTGTCGTCAG

MODE = seq_some *******************************************************************

>gi|152066|gb|M35122.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40987 seq_len:1040 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHAFIXA A.caulinodans nitrogen fixation protein (nifO and fixA) genes, complete cds and 5'end
>gi|152069|gb|M60872.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:40988 seq_len:1110 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 RHANODDA A.caulinodans nodD gene, complete cds
TCGTGCAGAGATACCATGCGCTGTGCGCCTACAGCGAAGCAAGATGTGACGCTGGACAATCTTTCGTAGCTACGGCAAAA
TTCACCGTGGGACGACACTGCGGATGCTGGGTAGATGGAGAGTTAGTTGCGATTTAAGGGACTTGATCTGAATCTGCTTG
CAGCGCAGAATATTCGGGAGCTTCCATATAGTCGTTTTCCGATATGCAGAGATCAAGAGCTCTCAATCGC
>gi|443814|emb|X77126.1|taxid:7 src_db:n kind: project: cat:B stage:nt src_type:o iid:262282 seq_len:982 divid:0 rank:species lineage:species=7,genus=6,family=335928,order=356,class=28211,phylum=1224,superkingdom=2,no_rank=131567,no_rank=1 AC16S6829 A.caulinodans (LMG 6829) 16S rRNA gene
"""


def test(out,mode):
    inp = StringIO(fastaInput)
    reader = FastaReader(inp)
    iRec = 0
    for rec in reader.records():
        taxid = int(rec.header().split("taxid:",1)[1].split(" ",1)[0])
        print "taxid = ", taxid
        out.write(rec.header())
        if mode == "seq":
            for line in rec.seqLines():
                out.write(line)
        elif mode == "seq_some":
            if iRec == 1:
                for line in rec.seqLines():
                    out.write(line)
        iRec += 1
                
    inp.close()

out = StringIO()
for mode in ("hdr","seq","seq_some"):
    out.write("\nMODE = "+mode+" *******************************************************************\n\n")
    test(out,mode)

outStr = out.getvalue()
print outStr
assert outStr == outExpected

out = StringIO()
test(out,"seq")
outStr = out.getvalue()
print outStr
assert outStr == fastaInput
