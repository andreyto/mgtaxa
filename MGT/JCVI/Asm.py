"""JCVI-specific sequence assembly processing"""

from MGT.Asm import *

from joblib import Memory
makedir("cache.tmp")
memory = Memory(cachedir="cache.tmp", verbose=0)

@memory.cache
def getSffIdReads(inp):
    return list(sffIdReadsReader(inp=inp))

@memory.cache
def getReadStatus(asmDir):
    return list(newblerReadStatusReader(asmDir=asmDir,format="row"))

def iterSffIdReads(iterSffTags):
    for (sff,tag) in iterSffTags:
        for id_read in getSffIdReads(inp=sff):
            yield (id_read,tag)

def jcviContigTagCount454(asmDir,tagExtractor,readContFilter=None):
    sffs = glob.glob(pjoin(asmDir,"sff","*.sff"))
    sffTags = [ (sff,tag) for (sff,tag) \
            in ( (sff,tagExtractor(sff)) for sff in sffs ) \
            if tag is not None ]
    idReadsDir = "id_reads"
    if False:
        makedir(idReadsDir)
        for sff,tag in sffTags:
            out = open(pjoin(idReadsDir,stripSfx(os.path.basename(sff))+".csv"),"w")
            for (id_read,tag) in iterSffIdReads(iterSffTags=((sff,tag),)):
                out.write("%s\n" % (id_read,))
            out.close()
    print "Searching these sff files and tags:"
    print "\n".join(("%s\t%s" % (sff,tag)) for (sff,tag) in sffTags)
    readContIter = newblerReadStatusReader(asmDir=asmDir,format="row")
    #readContIter = getReadStatus(asmDir)
    if readContFilter is not None:
        rcIter = readContFilter(readContIter)
    else:
        rcIter = readContIter
    return contigTagCount(readContIter=rcIter,
            readTagIter=iterSffIdReads(iterSffTags=sffTags))

