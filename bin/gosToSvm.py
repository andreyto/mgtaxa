from MGT.Common import *
from MGT.Svm import *
#from MGT.SampDb import *
#from MGT.SampDbKmer import *
#from MGT.PredictorDb import *



def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="append", type="string",dest="inSeq"),
        make_option("-o", "--out-dir",
        action="store", type="string",dest="outDir",default="."),
        make_option("-l", "--min-len",
        action="store", type="int",dest="minSampLen",default=750),
        make_option("-f", "--inp-format",
        action="store", type="choice",choices=("gos","ncbi","ca"),dest="inFormat",default="gos"),
        make_option("-e", "--degen-len",
        action="store", type="int",dest="degenLen",default=1),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

assert len(opt.inSeq) > 0, "Need at least one in-seq option"

def fastaToSvm(inFileFasta,outName,opt):
    svmWriter = SvmStringFeatureWriterTxt(outName+".samp")
    inpSeq = FastaReader(inFileFasta)
    if opt.degenLen >= 0:
        symCompr = SymbolRunsCompressor('N',opt.degenLen)
    else:
        symCompr = lambda s: s
    if opt.inFormat == "gos":
        meta, allLen = gosToSvm(inpSeq,svmWriter,symCompr,opt)
    elif opt.inFormat == "ca":
        meta, allLen = caToSvm(inpSeq,svmWriter,symCompr,opt)
    inpSeq.close()
    svmWriter.close()
    print "Saved %i samples out of %i total from file %s" % (len(meta.samp),len(allLen),inFileFasta)
    lenHist = numpy.histogram(allLen,bins=numpy.arange(0,allLen.max()+100,100,dtype='f8'))
    print "Original sample length histogram:\n%s\n%s" % lenHist
    dumpObj(meta,outName+".meta")

def gosToSvm(inpSeq,svmWriter,symCompr,opt):
    fldDtype=[("id","i8"),("id_mate","i8"),("forward",bool),("len_lib","i4"),("len_samp","i4")]
    fldValues = []
    allLen = []
    iRec = 0
    for rec in inpSeq.records():
        hdr = rec.header()
        parts = hdr.split()
        id_read = int(parts[0].split(">JCVI_READ_")[1])
        pairs = [ (pair[0].split('/')[1],pair[1]) for pair in 
                ( part.split('=') for part in parts[1:] ) ]
        pairs = dict(pairs)
        seq = symCompr(rec.sequence()[int(pairs["clr_range_begin"]):int(pairs["clr_range_end"])])
        lenSeq = len(seq)
        if lenSeq >= opt.minSampLen:
            svmWriter.write(0,seq,id_read)
            fldValues.append((id_read,
                        int(pairs["mate"]),
                        pairs["sequencing_direction"] == "forward",
                        int(pairs["library_id"].split("-")[-1].split("KB")[0])*1000,
                        lenSeq))
        allLen.append(lenSeq)
        if iRec % 1000 == 0:
            print "Processed %i records" % iRec
        iRec += 1
    sampMeta = numpy.asarray(fldValues,dtype=fldDtype)
    allLen = numpy.asarray(allLen,dtype='i4')
    meta = Struct(samp=sampMeta)
    return meta,allLen

def caToSvm(inpSeq,svmWriter,symCompr,opt):
    fldDtype=[("id","i8"),("len_samp","i4")]
    fldValues = []
    allLen = []
    iRec = 0
    for rec in inpSeq.records():
        hdr = rec.header()
        id = int(hdr.split(">scf")[1])
        seq = symCompr(rec.sequence())
        lenSeq = len(seq)
        if lenSeq >= opt.minSampLen:
            svmWriter.write(0,seq,id)
            fldValues.append((id,
                        lenSeq))
        allLen.append(lenSeq)
        if iRec % 1000 == 0:
            print "Processed %i records" % iRec
        iRec += 1
    sampMeta = numpy.asarray(fldValues,dtype=fldDtype)
    allLen = numpy.asarray(allLen,dtype='i4')
    meta = Struct(samp=sampMeta)
    return meta,allLen

#/usr/local/projects/GOSII/ANNOTATION/GSIOVIR110-I-01-3-4KB/GSIOVIR110-I-01-3-4KB.fasta
for inFile in opt.inSeq:
    outName = os.path.basename(inFile)
    outName = stripSfx(outName,'.')
    outName = os.path.join(opt.outDir,outName)
    assert not isSamePath(inFile,outName)
    makedir(outName)
    fastaToSvm(inFileFasta=inFile,outName=os.path.join(outName,"inp"),opt=opt)

