"""Utilities to extract information from sequence assemblies"""

from MGT.Common import *

def contigReadCount454(asmDir=None,inp=None,out=None):
    """Get per-contig read counts from 454ReadStatus.txt produced by Newbler assembler.
    We count 5' and 3' ends as 1/2 each to handle split reads.
    @param asmDir directory name for 454 assembly
    @param inp if not None, input stream in 454ReadStatus.txt format, otherwise it will be
    defined as asmDir/454ReadStatus.txt
    @param out output stream, if None - standard output will be used
    """
    inpClose = False
    if inp is None:
        inp = open(pjoin(asmDir,"454ReadStatus.txt"),'r')
        inpClose = True

    if out is None:
        out = sys.stdout

    cnt = defdict(float)

    # skip header
    inp.next()
    # count contigs, including cases when read spans two contigs
    # (we assign 1/2 to each - what else to do?)
    # FRDNTY201A47E7  Assembled       contig208724    402     -       contig134302    1772    +
    for line in inp:
        words = line.strip().split()
        if len(words) >= 3:
            cnt[words[2]] += 0.5
            cnt[words[5]] += 0.5
    if inpClose:
        inp.close()

    for (contig,count) in sorted(cnt.iteritems()):
        out.write("%s\t%i\n" % (contig,round(count)))

