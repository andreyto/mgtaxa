### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Entrez import *

#db=sequences
#http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&seq_start=1&seq_stop=9

req = EzRequest()

#for rec in req.fetchParsed([22855148,1068]):
for rec in req.fetchParsed(nrnd.randint(0,1000,100)):
    print rec
    for feat in rec.features:
        print feat.qualifiers.items()

