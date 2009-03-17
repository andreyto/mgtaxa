"""Application interface to graphical representation of real valued features."""
from MGT.FastaIO import *
from MGT.Svm import *
from MGT.App import *
from MGT.DirStore import *

__all__ = ["SeqImportApp"]


class SeqImportApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=("default",),
            dest="mode",default="default"),
            make_option("-i", "--in-seq",
            action="append", type="string",dest="inSeq"),
            make_option("-o", "--out-feat",
            action="store", type="string",dest="outFeat",default="samp"),
            make_option("-l", "--min-len",
            action="store", type="int",dest="minSampLen",default=0),
            make_option("-f", "--inp-format",
            action="store", type="choice",choices=("gos","ncbi","ca"),dest="inFormat",default="gos"),
            make_option("-e", "--degen-len",
            action="store", type="int",dest="degenLen",default=10),
            make_option("-p", "--prefix-id",
            action="store_true",dest="prefixId",default=True),
        ]
        return Struct(usage = "Import sequence samples for classification\n"+\
                "%prog [options]",option_list=option_list)

    def doWork(self,**kw):
        opt = self.opt
        print "App options:\n", opt
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        outFeat = self.store.getFilePath(opt.outFeat)
        for inFile in opt.inSeq:
            assert not isSamePath(inFile,outFeat)
        svmWriter = SvmStringFeatureWriterTxt(outFeat)
        for inFile in opt.inSeq:
            self.fastaToSvm(inFileFasta=inFile,svmWriter=svmWriter,opt=opt)
        svmWriter.close()


    def fastaToSvm(self,inFileFasta,svmWriter,opt):
        inpSeq = FastaReader(inFileFasta)
        if opt.degenLen >= 0:
            symCompr = SymbolRunsCompressor('N',opt.degenLen)
        else:
            symCompr = lambda s: s
        if opt.inFormat == "gos":
            meta, allLen = self.gosToSvm(inpSeq,svmWriter,symCompr,opt)
        elif opt.inFormat == "ca":
            meta, allLen = self.caToSvm(inpSeq,svmWriter,symCompr,opt)
        inpSeq.close()
        print "Saved %i samples out of %i total from file %s" % (len(meta.samp),len(allLen),inFileFasta)
        lenHist = numpy.histogram(allLen,bins=numpy.arange(0,allLen.max()+100,100,dtype='f8'))
        print "Original sample length histogram:\n%s\n%s" % lenHist
        ## @todo move saving of meta outside of loop or make file name unique
        #self.store.saveObj(meta,opt.outFeat+".meta")

    def gosToSvm(self,inpSeq,svmWriter,symCompr,opt):
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

    def caToSvm(self,inpSeq,svmWriter,symCompr,opt):
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


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(SeqImportApp)

