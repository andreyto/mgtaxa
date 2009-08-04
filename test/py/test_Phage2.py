### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Phage import *
from MGT.Taxa import *
taxaTree = loadTaxaTreeNew(allNames=True)
hosts = loadHosts(inpFile="viralHosts.pkl",taxaTree=taxaTree)
bacNode = taxaTree.getNode(bacTaxid)

if False:
    for node in taxaTree.iterDepthTop(virTaxid):
        name = node.name.lower()
        nameWords = name.split()
        if 'phage' in name and 'phage' not in nameWords:
            print "No phage word: ", node.name, node.names, node.id, node.rank, node.lineageStr()

for node in taxaTree.iterDepthTop(virTaxid):
    if hasattr(node,'hosts') and sum([h.node.isSubnode(bacNode) for h in node.hosts]): 
        if 'environmental samples' in node.name.lower():
            print "Env: ", node.id,node.name, node.rank, node.names,node.lineageStr(),'\n', node.hosts
        if len([ h for h in node.hosts if h.src == 'annot' ]) > 1:
            print "Annot: ", node.id,node.name, node.rank, node.names,node.lineageStr(),'\n', node. hosts
