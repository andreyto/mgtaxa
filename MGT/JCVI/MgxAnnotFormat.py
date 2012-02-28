"""Readers for output files of JCVI metagenomic annotation pipeline"""
from MGT.Common import *
import csv
from MGT.Bits.OrderedDict import *
from MGT.FastaIO import *
from MGT.Functional import *

_blastBtabDescr="""\
JCVI btab BLAST output format

1              Query Sequence Name

2              Date of the Analysis

3              Query Sequence Length

4              Search Method  --  Blast family application name

5              Database Name

6              Subject Sequence Name  --  Database entry name

7,8          Query Left End, Query Right End  --  The endpoints of the part
                of the query sequence which Blast aligns with the subject sequence.

9,10        Subject Left End, Subject Right End  --  The endpoints of the part
                of the subject sequence which Blast aligns with the query sequence.

11           Percent Identity  --  The fraction of residues which are absolute
                matches between the query and subject sequence, expressed in
                percent.  This number is calculated from the matching, similar,
                and aligned residues counted in the alignment.  It is checked
                against the stated values parsed out of the alignment header, but
                has greater precision.

12           Percent Similarity  --  The fraction of residues which are exact or
                similar matches between the query and subject sequence, expressed in
                percent.  This number is calculated from the matching and aligned
                residues counted in the alignment.  It is checked against the stated
                values parsed out of the alignment header, but has greater precision.
                It is only calculated for alignments which are based on amino acid
                sequences.

13           Score  --  The score assigned to the alignment by Blast

14,15     File Offset for Beginning of Alignment
                File Offset for End of Alignment -- Beginning and ending bytes of the
                portion of the Blast output file which actually shows the alignment.
                For Btab versions 2 and 3, this block contains the query sequence,
                alignment indicator characters, and subject sequence, but not the
                header information for an alignment.  Count is in bytes from
                beginning of Blast output file.  From version 4, a compile-time
                switch of MODIFIED_ALIGNMENTS is provided.  With this switch, the
                first pointer points to the start of the line containing the subject
                description. For multiple 'hits' for the same subject, last alignment
                for multiple hits will have the pointer to the description line for
                the first and the end of the alignment for the last hit.

16           Description  --  A freeform text field which contains the biological
                description field from the database for the subject sequence.  If
                this text occupies more than one line in the Blast output file, the
                NewLines are replaced by spaces.  Commas may occur in this field even
                if they are the field separator character, because this is the last
                field in the record.

17           Frame  --  1 through 6, or NULL

18           Query Strand  --  Plus, Minus or NULL

19           DB sequence length

20           Expect -- expected value - e value
"""

_blastBtabRecSample="""\
contig30487_725_2453_plus	Feb 10 2012	573	blastp	unspecified	UniRef100_D5CMJ6	1	369	2	370	65.04	80.76	1296.0	503.827		Filamentation induced by cAMP protein Fic n=1 Tax=Sideroxydans lithotrophicus ES-1 RepID=D5CMJ6_SIDLE	1	Plus	374	3.66518E-140	"""

def parseOrfOnContigId(id):
    ret = {}
    ret["id_cont"],ret["start_orf"],ret["end_orf"],ret["strand_orf"] = id.strip().split('_')
    ret["start_orf"] = int(ret["start_orf"]) - 1 # unit-offset is input, zero-offset is output
    ret["end_orf"] = int(ret["end_orf"])
    ret["strand_orf"] = 1 if ret["strand_orf"] == "plus" else -1
    return ret

class BtabReaderRecord(object):

    def __init__(self,rec):
        self.rec = rec
        self.parsed = rec.copy()

    def parse(self,fieldNames=None):
        rec = self.rec
        parsed = self.parsed
        if fieldNames is None:
            fieldNames = BtabReader.fieldNames
        map = BtabReader.formatMap
        for field in fieldNames:
            parsed[field] = map[field](rec[field])
        return parsed

class BtabReader(object):

    strC = lambda x: x.strip()
    pct = lambda x: float(x)/100.
    strandMapper = {"Plus":1,"Minus":-1,"NULL":0}

    formatMap = OrderedDict((
    ("id_q",strC),
    ("date",strC),
    ("len_q",int),
    ("method",strC),
    ("name_db",strC),
    ("id_s",strC),
    ("start_q",int),
    ("end_q",int),
    ("start_s",int),
    ("end_s",int),
    ("ident",pct),
    ("simil",pct),
    ("score",float),
    ("start_f_aln",strC), #this comes as a real number
    ("end_f_aln",strC), #this is empty
    ("descr",strC),
    ("frame",int),
    ("strand",lambda x,strandMapper=strandMapper: strandMapper[x]),
    ("len_s",int),
    ("e_val",float)
    ))

    fieldNames = formatMap.keys()
    
    def __init__(self,inp):
        x = openCsv(inp,
                mode="r",
                factory=csv.DictReader,
                fieldnames=self.fieldNames,
                dialect="excel-tab")
        self.reader = x.csvFile
        self.csvClose = x.csvClose
        self.csvFileInp = x.csvFileStream
        self.selRecs = {}

    def next(self):
        return BtabReaderRecord(self.reader.next())
    
    def __iter__(self):
        return self

    def close(self):
        if self.csvClose:
            self.csvFileInp.close()
           
    @coroutine
    def filterByIdQ(self,recProc=None):
        if recProc is None:
            recProc = lambda iter: iter
        idQ = (yield)
        for id_q,recs in it.groupby(lambda rec: rec["id_q"],self):
            if id_q == idQ:
                idQ = (yield recProc(recs))

class MappedPepFastaReader(FastaReader):
    """Reader for JCVI peptide FASTA file.
    Parses peptide-specific defline."""

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
        for rec in self.records():
            self.hdrParsed = rec.headerParsed()
            if self.hdrParsed["pep_id"] == idPep:
                return self.hdrParsed
        return None


_cameraAnnotRecSamp = """\
JCVI_PEP_metagenomic.orf.1132486264118.1	Hypothetical	N/A	1	common_name	hypothetical protein	gene_symbol		GO		EC		CAZY		TAXON	
JCVI_PEP_metagenomic.orf.1132486264118.1	PFAM::FullLength::Domain	PF03796	2	common_name	DnaB-like helicase C terminal domain	gene_symbol		GO	GO:0003824 || GO:0005524 || GO:0008152	EC		CAZY		TAXON	
JCVI_PEP_metagenomic.orf.1132486264118.1	PFAM::FullLength::Uncategorized	PF06745	1	common_name	KaiC	gene_symbol		GO		EC		CAZY		TAXON	
JCVI_PEP_metagenomic.orf.1132486264118.1	UnirefBLASTP::HighConfidence	UniRef100_C1HPD7	1	common_name	DNA repair protein radA	gene_symbol		GO	GO:0004252 || GO:0046872 || GO:0003684 || GO:0006508 || GO:0006281 || GO:0005524 || GO:0004176	EC		CAZY		TAXON	Escherichia
JCVI_PEP_metagenomic.orf.1132486264118.1	UnirefBLASTP::Reviewed	UniRef100_P24554	9	common_name	DNA repair protein radA	gene_symbol		GO	GO:0004252 || GO:0046872 || GO:0003684 || GO:0005524 || GO:0006508 || GO:0006281 || GO:0004176	EC		CAZY		TAXON	Escherichia coli
"""

class CameraAnnotRecord(object):

    #This is only used to assign partial ordering for selecting one best hit for each
    #ORF: the absolute value or difference are meaningless
    confMapper = {
            "Reviewed":6,
            "HighConfidence":5,
            "Putative":4,
            "ConservedDomain":3,
            "LowConfidence":2,
            "LowestConfidence":1,
            "Fragment":2
            }

    def __init__(self,rec):
        self.rec = rec
        self.parsed = {"id_q":rec[0].strip()}

    def recordType(self):
        parsed = self.parsed
        recType = self.rec[1].strip()
        try:
            recType,conf = self.rec[1].strip().split("::")
        except:
            parsed["type"] = recType
        else:
            if recType == "UnirefBLASTP": #prok pipeline
                recType = "UNIREF_BLASTP"
            parsed.update({
                "type" : recType,
                "conf" : conf,
                "conf_x" : self.confMapper[conf] if recType == "UNIREF_BLASTP" else 0.0
                })
        return recType

    def parse(self):
        parsed = self.parsed
        if not "type" in parsed:
            recType = self.recordType()
        else:
            recType = parsed["type"]
        rec = self.rec
        if recType == "UNIREF_BLASTP":
            parsed["id_s"] = rec[2].strip().split(":")[-1].strip() #prok pipeline has no prefix, vir does
            parsed["descr"] = rec[5].strip()
            parsed["tax_name"] = rec[15]
        return parsed

class CameraAnnotReader(object):

    def __init__(self,inp):
        x = openCsv(inp,
                mode="r",
                factory=csv.reader,
                dialect="excel-tab")
        self.reader = x.csvFile
        self.csvClose = x.csvClose
        self.csvFileInp = x.csvFileStream

    def next(self):
        return CameraAnnotRecord(self.reader.next())
    
    def __iter__(self):
        return self

    def close(self):
        if self.csvClose:
            self.csvFileInp.close()
           

class CameraAnnotReaderPepIds(CameraAnnotReader):
    """Specialization of @see CameraAnnotReader that maps peptide IDs to read IDs and coords"""

    orientationMapper = {-1 : "minus", 1 : "plus", 0 : "null"}

    def __init__(self,inp,pepMap):
        CameraAnnotReader.__init__(self,inp)
        self.pepMap = pepMap

    def next(self):
        ret = CameraAnnotRecord(self.reader.next())
        id_q = ret.parsed["id_q"]
        pepRec = self.pepMap.findIdPep(id_q)
        assert pepRec,"Peptide record not found in peptide map for peptide id: %s" % (id_q,)
        strand = self.orientationMapper[pepRec["orientation"]]
        ret.parsed["id_q"] = "%s_%s_%s_%s" % \
                (pepRec["read_id"],
                pepRec["begin"]+1,
                pepRec["end"],
                strand)
        ret.parsed["pep_id"] = pepRec["pep_id"]
        return ret
    
