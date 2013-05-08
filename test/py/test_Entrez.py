### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Entrez import *

#db=sequences
#http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&seq_start=1&seq_stop=9

req = EzRequest(rettype="genbank")

#for rec in req.fetchParsed([22855148,1068]):
for rec in req.fetch_blocked_parsed(nrnd.randint(0,1000,100)):
    print rec
    #for feat in rec.features:
    #    print feat.qualifiers.items()
    pass

ids = [
        33504458,
        33504458,
        363399092,
        307608751
        ]

req = EzRequest(rettype="xml")
for rec in req.summary(ids):
    print rec

