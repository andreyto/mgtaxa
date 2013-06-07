from MGT.Common import *
import csv
from MGT.Bits.OrderedDict import *
from MGT.FastaIO import *

def parseOrfOnContigId(id):
    ret = {}
    ret["id_cont"],ret["start_orf"],ret["end_orf"],ret["strand_orf"] = id.strip().split('_')
    ret["start_orf"] = int(ret["start_orf"]) - 1 # unit-offset is input, zero-offset is output
    ret["end_orf"] = int(ret["end_orf"])
    ret["strand_orf"] = 1 if ret["strand_orf"] == "plus" else -1
    return ret

class MappedPepFastaReader(FastaReader):
    """Reader for JCVI peptide FASTA file.
    Parses peptide-specific defline."""

    def __init__(self,*l,**kw):
        """Ctor.
        @param cacheParsed Flag to cache previously read parsed headers [False].
        Set to True if you expect to query records with findIdPep() in the order
        different from the order in the input file. Beware of memory consumption
        if you set this to True.
        """
        kw_c = kw.copy()
        cacheParsed = kw_c.pop("cacheParsed",False)
        FastaReader.__init__(self,*l,**kw_c)
        self.cacheParsed = cacheParsed
        if cacheParsed:
            self.parsedCache = {}
        else:
            self.parsedCache = None

    def headerParsed(self):
        """Return current header as a dict with parsed fields.
        Parses a single line like:
        >JCVI_PEP_metagenomic.orf.1132484412212.1 
        /orf_id=JCVI_ORF_metagenomic.orf.1132484412211.1 
        /read_id=deg120007911192 
        /begin=0 /end=694 /orientation=-1 
        /5_prime_stop=0 /5_prime_loc=NA 
        /3_prime_stop=0 /ttable=11 
        /length=694 /read_defline="""
        hdr = self.header()
        hdr = hdr.split(">",1)[1]
        pep_id,rest = hdr.split(None,1)
        ret = dict(pep_id=pep_id.strip())
        #this will not handle another /something_defline="" inside
        #/read_defline="", but hopefully we will not see such input
        for tag,value in ( match.groups() for match \
                in re.finditer(r'(?:/(\S*)=((?:\"[^\"]*\")|(?:\S*\b)))+',rest) ):
            tag = tag.strip()
            assert tag not in ret,"Repeating tags in defline: %s" % (rest,)
            ret[tag] = value.strip().strip('"').strip()
        #coords seem to be between the bases or zero-offset
        for tag in ("begin","end","orientation","length"):
            if tag in ret:
                ret[tag] = int(ret[tag])
        return ret 
    
    def findIdPep(self,idPep):
        
        if hasattr(self,"hdrParsed") and self.hdrParsed["pep_id"] == idPep:
            return self.hdrParsed
        
        parsedCache = self.parsedCache

        if parsedCache is not None:
            parsedFound = parsedCache.get(idPep,None)
            if parsedFound:
                return parsedFound
                
        for rec in self.records():
            hdrParsed = rec.headerParsed()
            self.hdrParsed = hdrParsed
            if parsedCache is not None:
                parsedCache[hdrParsed["pep_id"]] = hdrParsed
            if hdrParsed["pep_id"] == idPep:
                return self.hdrParsed
        
        return None

