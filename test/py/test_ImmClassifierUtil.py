from MGT.ImmClassifierAppUtils import *

headers = """\
>gi|374264089|ref|NC_016756.1| Ashbya gossypii FDAG1 mitochondrion, complete genome
>gi|113200973|gb|DQ899947.1| Liriodendron tulipifera chloroplast, complete genome
>gi|113200973|gb|DQ899947.1| Liriodendron tulipifera something
>gi|374333421|ref|NC_016646.1| Pseudovibrio sp. FO-BEG1 plasmid unnamed, complete sequence"""


hdrParser = FastaHeaderParser()

for hdr in headers.splitlines():
    hdrNew = hdrParser.parseAndTagPristineHeader(hdr)
    print hdrNew
    res = hdrParser.parseTaggedHeader(hdrNew)
    print str(res)
    print GenomicElementType.typeName[res.genElType]

