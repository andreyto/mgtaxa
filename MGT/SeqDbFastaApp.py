### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application interface to SeqDbFasta class."""

__all__ = ["SeqDbFastaApp"]

from MGT.App import *
from MGT.Common import *
from MGT.SeqDbFasta import *
from MGT.FastaIO import *
from MGT.FastaSplitters import *

class SeqDbFastaApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", 
            type="choice",
            choices=("split",),
            dest="mode",
            help="split - create new store split into equal chunks"),
            
            optParseMakeOption_Path(None, "--inp-seq",
            dest="inpSeq",
            help="Path of the input store"),

            optParseMakeOption_Path(None, "--inp-seq-db-ids",
            dest="inpSeqDbIds",
            help="Optional JSON file with a selection of SeqDb IDs if inpSeq is a store"),

            optParseMakeOption_Path(None, "--out-seq-db-ids",
            dest="outSeqDbIds",
            help="Optional JSON file name to write all SeqDb IDs in output outSeq"),

            optParseMakeOption_Path(None, "--out-seq",
            dest="outSeq",
            help="Path of the output store"),
            
            make_option(None, "--chunk-size",
            action="store", 
            type="int",
            dest="chunkSize",
            help="Size of chunk in output store for --mode split"),

            make_option(None, "--seq-filter",
            action="store", 
            type="choice",
            choices=("len_degen","none"),
            dest="seqFilter",
            default="len_degen",
            help="Filter input with this filter type"),
            
            make_option(None, "--min-len-seq",
            action="store", 
            type="int",
            dest="minLenSeq",
            help="Minimum sequence length if len_degen filter is used"),
            
            make_option(None, "--assert-out-seq",
            action="store", 
            type="int",
            dest="assertOutSeq",
            default=1,
            help="Assert that output store is not empty"),
        ]
        return Struct(usage = "Operations on SeqDbFasta store.\n"+\
                "%prog [options]",option_list=option_list)
        
    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        pass
    
    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "split":
            return self.split(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def split(self,**kw):
        """Convert input sequence into split store, optionally filtering.
        Parameters are taken from self.opt
        """
        opt = self.opt
        to_close = []
        try:
            if SeqDbFasta.isStore(opt.inpSeq):
                inpSeq = SeqDbFasta.open(opt.inpSeq,"r")
                to_close.append(inpSeq)
                if opt.isDef("inpSeqDbIds"):
                    inpSeqDbIds = load_config_json(opt.inpSeqDbIds)
                else:
                    inpSeqDbIds = inpSeq.iterIds()
                reader = inpSeq.fastaReaderChain(inpSeqDbIds)
                to_close.append(reader)
            else:
                reader = FastaReaderChain( (inp for inp in glob.glob(opt.inpSeq)) )
                to_close.append(reader)
            with closing(SeqDbFasta.open(opt.outSeq,"c")) as outSeq:
                if opt.seqFilter is None or opt.seqFilter == "none":
                    nRec = splitFastaReaderIntoChunks(
                            reader=reader,
                            outStore=outSeq,
                            maxChunkSize=opt.chunkSize,
                            compresslevel=1)
                elif opt.seqFilter == "len_degen":
                    nRec = splitFastaReaderIntoChunksLengthDegen(
                            reader=reader,
                            outStore=outSeq,
                            maxChunkSize=opt.chunkSize,
                            minSeqLen=opt.minLenSeq,
                            compresslevel=1,
                            minNonDegenRatio=0.50)
                else:
                    raise ValueError("Unknown seqFilter value: {}".format(opt.seqFilter))
                if opt.isDef("assertOutSeq") or opt.isDef("outSeqDbIds"):
                    outSeqDbIds = outSeq.getIdList()
                if opt.assertOutSeq:
                    assert outSeqDbIds, \
                            ("Output sequence store {} has no "+\
                            "partitions after creating from: {}")\
                            .format(opt.outSeq,opt.inpSeq)
                    for (idSeqDb,lengths) in  outSeq.seqLengthsMany(outSeqDbIds):
                        if len(lengths) > 0:
                            break
                    else:
                        raise AssertionError(("No sequences remain in output store {}"+\
                                "created from: {}").format(opt.outSeq,opt.inpSeq))
                if opt.isDef("outSeqDbIds"):
                    save_config_json(outSeqDbIds,opt.outSeqDbIds)
                outSeq.setAttr("seq_count",nRec)
                outSeq.save()
        finally:
            #if one close() raises, remaining are not called
            while to_close:
                to_close.pop().close()

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(SeqDbFastaApp)

