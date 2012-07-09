"""Utilities to extract information from sequence assemblies"""

from MGT.Common import *
import collections

NewblerReadStatusRecordCol = collections.namedtuple('NewblerReadStatusRecordCol', 
        ["id_read", "id_cont1", "id_cont2"])

NewblerReadStatusRecordRow = collections.namedtuple('NewblerReadStatusRecordRow', 
        ["id_read", "id_cont", "i_cont"])

def newblerReadStatusReader(asmDir=None,inp=None,format="col"):
    """Get per-contig read counts from 454ReadStatus.txt produced by Newbler assembler.
    We count 5' and 3' ends as 1/2 each to handle split reads.
    @param asmDir directory name for 454 assembly
    @param inp if not None, input stream in 454ReadStatus.txt format, otherwise it will be
    defined as asmDir/454ReadStatus.txt
    @param format ["row"|"col"] - if row - return two rows per record, otherwise
    return one record with two contig columns, second row is returned only if
    contig is different
    """
    inpClose = False
    if inp is None:
        inp = open(pjoin(asmDir,"454ReadStatus.txt"),'r')
        inpClose = True

    # skip header
    inp.next()
    # count contigs, including cases when read spans two contigs
    # (we assign 1/2 to each - what else to do?)
    # FRDNTY201A47E7  Assembled       contig208724    402     -       contig134302    1772    +
    for line in inp:
        words = line.strip().split()
        if len(words) >= 3:
            if format == "col":
                yield NewblerReadStatusRecordCol(id_read=words[0],
                        id_cont1=words[2],
                        id_cont2=words[5])
            elif format == "row":
                yield NewblerReadStatusRecordRow(id_read=words[0],
                        id_cont=words[2],
                        i_cont=0)
                if words[2] != words[5]:
                    yield NewblerReadStatusRecordRow(id_read=words[0],
                            id_cont=words[5],
                            i_cont=1)
            else:
                raise ValueError("Unknown format parameter: %s" % (format,))

    if inpClose:
        inp.close()

def contigReadCount454(asmDir=None,inp=None,out=None):
    """Get per-contig read counts from 454ReadStatus.txt produced by Newbler assembler.
    We count 5' and 3' ends as 1/2 each to handle split reads.
    @param asmDir directory name for 454 assembly
    @param inp if not None, input stream in 454ReadStatus.txt format, otherwise it will be
    defined as asmDir/454ReadStatus.txt
    @param out output stream, if None - standard output will be used
    """
    if out is None:
        out = sys.stdout

    cnt = defdict(float)
    
    for rec in newblerReadStatusReader(asmDir=asmDir,inp=inp):
        cnt[rec.id_cont1] += 0.5
        cnt[rec.id_cont2] += 0.5

    for (contig,count) in sorted(cnt.iteritems()):
        out.write("%s\t%i\n" % (contig,round(count)))

def sffIdReadsReader(inp):
    cmd = ["sffinfo","-a",inp]
    if not os.path.isfile(inp):
        raise OSError("Input file does not exist: %s", inp)
    p = Popen(cmd,env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True)
    if p.returncode:
        raise CalledProcessError(str(cmd),p.returncode)
    for line in p.stdout:
        yield line.strip()
    p.stdout.close()
    p.wait()

def contigTagCount(readContIter,readTagIter):
    """Get per-contig tag counts.
    @param readCont iterator of (id_read,id_cont) tuples (can be N-to-M relation as 
    in Newbler output where one read can be in two contigs)
    @param readTagIter input iterator of (id_read,tag) tuples
    @return dict(id_contig->dict(tag->count))

    Example use case is when tags are library IDs, and we want to get
    library composition of every contig.
    """
    contCache = ObjCache()
    #readToCont = defdict(set)
    readToCont = dict() 
    for rec in readContIter:
        id_read = rec[0]
        id_cont = contCache(rec[1])
        if not id_read in readToCont:
            readToCont[id_read] = id_cont
        else:
            if hasattr(readToCont[id_read],"add"):
                readToCont[id_read].add(id_cont)
            else:
                readToCont[id_read] = set((readToCont[id_read],id_cont))

    #readToCont = dict(readToCont)
    #the line below fails on a mix of sets and strings in values
    #print "\n".join(("DEBUG:\tRC\t%s\t%s" % (x,y)) \
    #        for (x,y) in sorted(readToCont.items(),key=lambda v: (v[1],v[0])))

    contTagCnt = defdict(lambda: defdict(int))
    
    for (id_read,tag) in readTagIter:
        if id_read in readToCont:
            cont = readToCont[id_read]
            if hasattr(cont,"add"):
                for id_cont in cont:
                    contTagCnt[id_cont][tag] += 1
                    print "DEBUG:\tRCT\t%s\t%s\t%s" % (id_read,id_cont,tag)
            else:
                contTagCnt[cont][tag] += 1
                print "DEBUG:\tRCT\t%s\t%s\t%s" % (id_read,cont,tag)
    contCache.clear()
    contTagCnt = dict(contTagCnt)
    for key in contTagCnt:
        contTagCnt[key] = dict(contTagCnt[key])

    return contTagCnt

def contigTagCount454(asmDir=None,inp=None,readTagIter=None):
    """Get per-contig tag counts from 454ReadStatus.txt produced by Newbler assembler.
    @param asmDir directory name for 454 assembly
    @param inp if not None, input stream in 454ReadStatus.txt format, otherwise it will be
    defined as asmDir/454ReadStatus.txt
    @param readTagIter input iterator of (id_read,tag) tuples
    @return dict(id_contig->dict(tag->count))

    Example use case is when tags are library IDs, and we want to get
    library composition of every contig.
    """
    readContIter = newblerReadStatusReader(asmDir=asmDir,inp=inp,format="row")
    return contigTagCount(readContIter=readContIter,readTagIter=readTagIter)

def contigReadCountVelvet(contFasta,kmerSize,readLen,out=None):
    """Get per-contig read counts from Velvet contig deflines.
    Velvet manual describes the defline format as well as the
    relation between k-mer coverage provided in the defline and
    base coverage.
    @param contFasta path to Velvet contig file (usually called contigs.fa).
    If None, standard input is read
    @param kmerSize k-mer size (required to compute base coverage)
    @param readLen
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

