### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Svm import *
from MGT.FastaIO import FastaReader
#from MGT.SampDb import *
#from MGT.SampDbKmer import *
#from MGT.PredictorDb import *



def getProgOptions(args):
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="store", type="string",dest="inSeq"),
        make_option("-o", "--out-seq",
        action="store", type="string",dest="outSeq"),
        make_option("-l", "--min-len",
        action="store", type="int",dest="minSampLen",default=0),
        make_option("-f", "--inp-format",
        action="store", type="choice",choices=("gos","ncbi","ca","fasta"),dest="inFormat",default="fasta"),
        make_option("-z", "--out-format",
        action="store", type="choice",choices=("svm","fasta"),dest="outFormat",default="fasta"),
        make_option("-e", "--degen-len",
        action="store", type="int",dest="degenLen",default=-1),
        make_option("-p", "--prefix-id",
        action="store_true",dest="prefixId",default=False),
        make_option(None, "--fasta-line-len",
        action="store", type="int",dest="fastaLineLen",default=1000),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args(args=args)

    return options,args


def fastaToSvm(inFileFasta,outName,opt):
    assert not isSamePath(inFileFasta,outName)
    if opt.outFormat == "svm":
        svmWriter = SvmStringFeatureWriterTxt(outName)
    elif opt.outFormat == "fasta":
        svmWriter = SvmFastaFeatureWriterTxt(outName,lineLen=opt.fastaLineLen)
    inpSeq = FastaReader(inFileFasta)
    if opt.degenLen >= 0:
        symCompr = SymbolRunsCompressor('N',opt.degenLen)
    else:
        symCompr = lambda s: s
    if opt.inFormat == "gos":
        meta, allLen = gosToSvm(inpSeq,svmWriter,symCompr,opt)
    elif opt.inFormat == "ca":
        meta, allLen = caToSvm(inpSeq,svmWriter,symCompr,opt)
    else:
        meta, allLen = genericFastaToSvm(inpSeq,svmWriter,symCompr,opt)
    inpSeq.close()
    svmWriter.close()
    print "Saved %i samples out of %i total from file %s" % (len(meta.samp),len(allLen),inFileFasta)
    lenHist = numpy.histogram(allLen,bins=numpy.arange(0,allLen.max()+100,100,dtype='f8'))
    print "Original sample length histogram:\n%s\n%s" % lenHist
    dumpObj(meta,outName+".meta")

def gosToSvm(inpSeq,svmWriter,symCompr,opt):
    fldDtype=[("id","i8"),("id_mate","i8"),("forward",bool),("len_lib","i4"),("len_samp","i4"),("library_id","S40")]
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
            library_id = pairs["library_id"].split("JCVI_LIB_")[1]
            metaRec = ( id_read,
                        int(pairs["mate"]),
                        pairs["sequencing_direction"] == "forward",
                        int(pairs["library_id"].split("-")[-1].split("KB")[0])*1000,
                        lenSeq,
                        library_id )
            fldValues.append(metaRec)
            if opt.prefixId:
                id_read = "%s_%s" % (library_id,id_read)
            svmWriter.write(0,seq,id_read)
        allLen.append(lenSeq)
        if iRec % 1000 == 0:
            print "Processed %i records" % iRec
        iRec += 1
    sampMeta = numpy.asarray(fldValues,dtype=fldDtype)
    allLen = numpy.asarray(allLen,dtype='i4')
    meta = Struct(samp=sampMeta)
    return meta,allLen

def genericFastaToSvm(inpSeq,svmWriter,symCompr,opt):
    fldDtype=[("id",idDtype),("len_samp","i4")]
    fldValues = []
    allLen = []
    iRec = 0
    for rec in inpSeq.records():
        hdr = rec.header()
        id = rec.getId()
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

caToSvm = genericFastaToSvm

if __name__ == "__main__":
    opt,args = getProgOptions(args=None)
    assert opt.inSeq is not None, "Need at least --in-seq option"
    fastaToSvm(inFileFasta=opt.inSeq,outName=opt.outSeq,opt=opt)

