### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Phage import *

from MGT.Taxa import *

taxaTree = loadTaxaTreeNew(allNames=True)
phParser = VirHostParser(taxaTree=taxaTree)
#phParser.assignHostsByVirName()
#phParser.assignHostsByGbFeat(gbFile="refseq/viral.genomic.gbff")
phParser.assignHostsAndSave(gbFile="refseq/viral.genomic.gbff.hdr",outFile="viralHosts.pkl")
hosts = loadHosts(inpFile="viralHosts.pkl",taxaTree=taxaTree)
#db=sequences
#http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=5&rettype=gb&seq_start=1&seq_stop=9
