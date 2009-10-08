"""Support for writing GFF3 files"""

class GFF3Header(object):
    """GFF3 file header"""
    
    def __str__(self):
        return "##gff-version 3\n"


class GFF3Attributes(dict):
    """Represents GFF3 attributes."""
    def __str__(self):
        return ';'.join(( "%s=%s" % (tag, value if isinstance(value,str) or \
                not hasattr(value,"__len__") else ','.join( ("%s" % x for x in value) )) for \
                (tag,value) in sorted(self.items())))

GFF3At = GFF3Attributes

class GFF3Record(object):
    """Represents one record (line) in GFF3 file."""


    def __init__(self,seqid='.',source='.',type='.',start=0,end=1,
            score='.',strand='.',phase='.',attribs=None,**kw):
        """Constructor.
        Attributes can be passed as GFF3Attributes object through attribs argument, or
        as additional keyword arguments, or both.
        @param start Start of feature (zero-based, will be converted to GFF unit-based on output)
        @end end of feature (zero-based, will be converted to GFF unit-based on output)"""
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        if attribs is None:
            self.attribs = GFF3Attributes(kw)
        else:
            self.attribs = GFF3Attributes(attribs)
            self.attribs.update(kw)

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.seqid,self.source,self.type,self.start+1,self.end,
                self.score,self.strand,self.phase,self.attribs)

    def copy(self):
        from copy import copy
        c = copy(self)
        c.attribs = copy(self.attribs)
        return c

    def fromSeqFeat(self,feat):
        """Pull data from Bio.SeqFeature instance"""
        self.type = feat.type
        self.start = feat.location.nofuzzy_start
        self.end = feat.location.nofuzzy_end
        self.strand = '+' if feat.strand > 0 else '-' if feat.strand < 0 else '.'
        ats = self.attribs
        quals = feat.qualifiers
        ret = [ self ]
        if "id" in quals:
            ats["ID"] = quals["id"][0]
        if self.type == "CDS":
            self.phase = 0
            ats["Name"] = quals["product"][0]
            ats["ID"] = quals["protein_id"]
        elif self.type == "gene":
            ats["Name"] = quals["locus_tag"][0]
            ats["ID"] = ats["Name"]
        else:
            ats["Name"] = quals.get("note",quals.get("id",(self.type,)))[0]
        if quals.get("rpt_type",("",))[0] == "CRISPR":
            self.type = "CRISPR"
            self.strand = '.'
            ats["ID"] = quals["note"]
            if "rpt_unit_range" in quals:
                for srange in quals["rpt_unit_range"]:
                    start,end = srange.split("..")
                    start = int(start)-1
                    end = int(end)
                    rec = self.copy()
                    rec.type = "repeat_unit"
                    rec.start = start
                    rec.end = end
                    rec_ats = rec.attribs
                    del rec_ats["Name"], rec_ats["ID"]
                    rec_ats["Parent"] = ats["ID"]
                    ret.append(rec)
        return ret
                    


