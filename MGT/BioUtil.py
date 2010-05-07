
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Extensions for BioPython.
"""


from MGT.Common import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation

def extractFeatureLocIndex(feat,format="range"):
    """Return indices of every position that belongs to SeqFeature instance.
    It uses nofuzzy_end and nofuzzy_start approximations. This might or might not be an issue
    for applications like extracting non protein coding nucleotide sequences.
    @param feat SeqFeature instance
    @param format return format
    @return tuple of two arrays of equal length: first is array of positions, second if array of strand values
    if 'format' is "range" [default], positions are expressed as tuple ranges, if format is "ind",
    positions are indices of each base that belongs to this feature
    """
    assert format in ("range","ind"),"Unknown format argument: %s" % (format,)
    pos = []
    strand = []
    # converts strand None to 0
    _coalesce = lambda x: x if x else 0
    def _extract(self,par_strand):
        # The way GenBank parser (v1.53) works, recursive complement(...) operator is limited
        # Specifically, the inner complment() in complement(join(207497..208369,complement(1..687)))
        # will be silently ignored by the parser.
        self_strand = _coalesce(self.strand)
        assert par_strand >= 0 or self_strand == par_strand,"Mismatch in sub-feature strand:\n %s" % (feat,)
        if self.sub_features:
            if self.location_operator!="join":
                raise ValueError(self.location_operator)
            for f_sub in self.sub_features:
                _extract(f_sub,self_strand)
        else:
            pos.append((self.location.nofuzzy_start,self.location.nofuzzy_end))
            strand.append(self_strand)
    _extract(feat,_coalesce(feat.strand))
    if format == "range":
        return (n.asarray(pos,dtype=int),n.asarray(strand,dtype="i1"))
    elif format == "ind":
        return (n.concatenate([n.arange(p[0],p[1]) for p in pos]),
                n.concatenate([n.ones(p[0][1]-p[0][0],dtype="i1")*p[1] for p in it.izip(pos,strand)]))

class FeatureLocIndex(object):
    """Accumulator for feature locations"""

    def __init__(self):
        self.pos = []
        self.strand = []

    def extract(self,feat):
        """Extract location of a single feature.
        This is supposed to be called repeatedly, and the results are
        accumulated inside this instance."""
        pos,strand = extractFeatureLocIndex(feat,format="range")
        self.pos.append(pos)
        self.strand.append(strand)

    def toFasta(self,rec,outFeat=None,outNonFeat=None):
        """Output feature sequences and/or non-feature (e.g. non-coding regions) sequences as FASTA.
        This should be called after feature locations have been accumulated with calls to extract().
        @param rec SeqRecord on which the accumulated features are defined.
        @param idPrefix IDs in FASTA deflines will be constructed as "idPrefix_i_j" where i is feature
        index starting from 0 in the list of accumulated features, and j is location range starting from
        0 within a single feature
        @param outFeat Output stream for feature sequences
        @param outNonFeat Output stream for non-feature sequences.
        @note Features will be written with strand taken into account, but inter-feature regions will
        be writen always on the positive strand.
        @note The current implementation will merge overlapping features, with strand of the
        overlapping region picked from one of the features in undefined order"""
        from itertools import izip
        mask = n.zeros(len(rec.seq),dtype=bool)
        for (iFeat,(pos,strand)) in enumerate(izip(self.pos,self.strand)):
            for (iRange,(p,s)) in enumerate(izip(pos,strand)):
                mask[p[0]:p[1]] = True
                if outFeat:
                    seq = rec.seq[p[0]:p[1]]
                    if s < 0:
                        seq = seq.reverse_complement()
                    r = SeqRecord(id="%s_%s_%s" % (rec.id,iFeat,iRange),
                            description="Feature=%s Range=%s Strand=%s" % (iFeat,iRange,s),
                            seq=seq)
                    SeqIO.write([r],outFeat,"fasta")
        
        if outNonFeat:
            mask = n.logical_not(mask)
            pos = runsOfOnesArray(mask)
            iFeat = 0
            s = 0
            for iRange,p in enumerate(pos):
                seq = rec.seq[p[0]:p[1]]
                r = SeqRecord(id="%s_%s_%s" % (rec.id,iFeat,iRange),
                        description="Feature=%s Range=%s Strand=%s" % (iFeat,iRange,s),
                        seq=seq)
                SeqIO.write([r],outNonFeat,"fasta")

def extractCdsAndNonCds(inpGb,outCds=None,outNonCds=None):
    """Split sequences from nucleotide GenBank flat file into CDS and non-CDS FASTA streams.
    @param[in] inpGb GenBank file
    @param[out] Output stream for CDS (can be None to skip writing CDS)
    @param[out] Output stream for non-CDS (can be None to skipt writing non-CDS)
    """
    for rec in SeqIO.parse(inpGb,"genbank"):
        flind = FeatureLocIndex()
        for feat in rec.features:
            if feat.type == "CDS":
                flind.extract(feat)
        flind.toFasta(rec=rec,outFeat=outCds,outNonFeat=outNonCds)

