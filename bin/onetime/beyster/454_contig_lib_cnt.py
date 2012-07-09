"""Get per-contig site counts from 454ReadStatus.txt produced by Newbler assembler.
"""

from MGT.JCVI.Asm import *

def tagExtractorSite(sff):
    import re
    #extract frist 3 digits ('679') from file name e.g.
    #GOSEUROPE679-RLB-RL6-01-819_GLDFQNX01_RL006_ACGCGTCTAGT.trimmed.sff
    siteRx = r"^[A-Z]*([0-9]{3})"
    name = os.path.basename(sff).strip()
    tag = re.match(siteRx,name).group(1)
    return tag

def tagExtractorLib(sff):
    tagSite = tagExtractorSite(sff)
    tag = os.path.basename(sff).strip().split("_")[0].split("-")[0]
    return tag
    if tagSite in ("678","679","681"):
        return tag
    else:
        return None

def tagExtractorFile(sff):
    tagSite = tagExtractorSite(sff)
    tag = os.path.basename(sff).strip()
    return tag
    if tagSite in ("678","679","681"):
        return tag
    else:
        return None

asmDir = sys.argv[1]
if len(sys.argv) >= 3:
    id_contSet = set((x.strip() for x in open(sys.argv[2],'r').read().split()))
    def readContFilter(iterReadCont,id_contSet=id_contSet):
        for rec in iterReadCont:
            if rec[1] in id_contSet:
                yield rec
    rcFilt = readContFilter
else:
    rcFilt = None

contTagCnt = jcviContigTagCount454(asmDir=asmDir,
        tagExtractor=tagExtractorLib,
        readContFilter=rcFilt)

out = sys.stdout

out.write("id_cont\ttag\tcnt\n")

for id_cont in sorted(contTagCnt.keys()):
    for tag,cnt in sorted(contTagCnt[id_cont].items(),
            key=lambda v: (v[0],-v[1])):
        out.write("%s\t%s\t%s\n" % (id_cont,tag,cnt))





