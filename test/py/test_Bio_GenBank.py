### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from MGT.Util import *
from Bio import SeqIO
from Bio.Seq import Seq,UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

#gbFile = pjoin(options.testDataDir,"iofilt.test.gz")
gbFile = "/usr/local/projects/GOSII/syooseph/MF150/all_Moore_data/CP000878.1.gbk.gz"
inGb = ioFilter(openCompressed(gbFile,'r'),code="lambda x: x.replace('&gt;','').replace('&lt;','')",mode="line")
outGb = openCompressed("test.out.gb.gz",'w')
for rec in SeqIO.parse(inGb,"genbank"):
    ft = SeqFeature(location=FeatureLocation(200,300),
            type="repeat_region",
            strand=0,
            qualifiers={"rpt_type":("CRISPR",)},
            id="CRISPR01")
    rec.features.append(ft)
    #rec.seq = UnknownSeq(len(rec.seq),alphabet=IUPACAmbiguousDNA())
    SeqIO.write([rec],outGb,"genbank")
inGb.close()
outGb.close()

