"""Helper classes and methods for ImmClassifierApp"""
from MGT.Common import *

class GenomicElementType(object):

    typeEnum = {
            "other" : 0x00,
            "chromosome" : 0x01,
            "plasmid" : 0x02,
            "chloroplast" : 0x04,
            "rhodoplast" : 0x08,
            "mitochondrion" : 0x10,
            "virus" : 0x20,
            "phage" : 0x40
            }


    typeRe = (
            {"name":"chromosome"},
            {"name":"plasmid"},
            {"name":"chloroplast"},
            {"name":"rhodoplast"},
            {"name":"mitochondrion"},
            {"name":"virus","re":r"\S*virus"},
            {"name":"phage","re":r"\S*phage"},
            )

    for el in typeRe:
        el["re"] = re.compile(r"\b"+el.get("re",el["name"])+r"\b",re.I)

    _gel = typeEnum
    _gel["plastid"] = _gel["chloroplast"] \
            | _gel["rhodoplast"] \
            | _gel["mitochondrion"]
    _gel["extra_chrom"] = _gel["plasmid"] \
            | _gel["plastid"]

    _gel["viral"] = _gel["virus"] \
            | _gel["phage"]
    
    typeName = dict((enum,name) for (name,enum) in typeEnum.items())


class FastaHeaderParser(object):
    gElType = GenomicElementType()
    tag_prefix = "mgt_tag_"
    gelt_tag_name = "gelt"

    def parsePristineHeader(self,hdr):
        rexps = self.gElType.typeRe
        enums = self.gElType.typeEnum
        enRes = enums["other"]
        for rexp in rexps:
            if rexp["re"].search(hdr):
                enRes |= enums[rexp["name"]]
        return Struct(genElType=enRes)
    
    def tagPristineHeader(self,hdr,parsed):
        tag = self.tag_prefix+self.gelt_tag_name
        assert tag not in hdr,"Found my tag in supposedely untouched FASTA header"
        hdr = "%s %s=%#x\n" % (hdr.rstrip(),tag,parsed.genElType)
        return hdr

    def parseAndTagPristineHeader(self,hdr):
        parsed = self.parsePristineHeader(hdr)
        return self.tagPristineHeader(hdr,parsed)

    def parseTaggedHeader(self,hdr):
        tag = self.tag_prefix+self.gelt_tag_name+"="
        val = int(hdr.split(tag)[1].split()[0],16)
        return Struct(genElType=val)

class FastaTrainingSeqFilter(object):
    """Filter for FastaReader that implements policy for selecting genomic type elements"""

    def __init__(self,policy):
        self.policy = policy
        self.hdrParser = FastaHeaderParser()

    def __call__(self,hdr,seq):
        hdrParser = self.hdrParser
        typeEnum = hdrParser.gElType.typeEnum
        parsed = hdrParser.parsePristineHeader(hdr)
        action = "pass"
        if self.policy == "drop-plasmids":
            if parsed.genElType & typeEnum["plasmid"]:
                action = "drop"
        elif self.policy == "extra-chrom-only":
            if not parsed.genElType & typeEnum["extra_chrom"]:
                action = "drop"
        else:
            raise ValueError("Unknown value for %s" % (self.policy,))
        if action == "pass":
            return hdrParser.tagPristineHeader(hdr,parsed),seq
        else:
            return None

