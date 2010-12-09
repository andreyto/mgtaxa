### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Classes and methods to build a database of known pairs of bacteriophages and their microbial hosts. 
This is done by first parsing the available 'host' assignments from NCBI GenBank records, and for the remaining 
unassigned phages - by extracting the microbial names from the phage names."""

__all__ = ["VirHostParser","PhageHostSeqPicker","loadHosts"]

from MGT.TaxaTree import *
from Bio import SeqIO
from Svm import IdLabels

def splitNameToHosts(name):
    """Take a full name, and return a list of diminishing prefixes before 'phage' or 'virus'.
    @param name name of virus or host, e.g.'blackberry cv. apache virus'
    @ret list of (num_of_removed_parts,prefix) pairs, e.g. [(0,'blackberry cv. apache'), (1,'blackberry cv.'), (2,'blackberry')]"""
    name = name.lower()
    words = name.split()
    hostWords = []
    if words[0] in ("uncultured","unclassified","unidentified","hybrid"):
        words = words[1:]
    for word in words:
        if "phage" in word:
            # this is an important kind for us, so do it explicitely
            if 'cyanophage' == word:
                hostWords.append('cyanobacteria') #phylum
            else:
                stem = word.split("phage")[0]
                if stem == 'bacterio': #word was 'bacteriophage'
                    pass
                elif stem == 'viro': #word was 'virophage'
                    pass
                elif len(stem) > 0:
                    #'' - as in 'vibriophage', 'bacter' and 'bacterium' - genus, 'bacteria' - phylum
                    # we can go through subtracting/adding the full list of official rank suffixes
                    # http://www.bacterio.cict.fr/classification.html
                    # but the ones used below seem to cover all current records
                    # Also, the suffix rules cover only ranks above genus and up to class inclusive,
                    # and use genus name as a stem. Above classes, there are no rules.
                    hostWords.append([ stem+sfx for sfx in ('','bacter','bacterium','bacteria') ])
                else:
                    pass
            break
        elif "virus" in word:
            break
        else:
            hostWords.append(word)
    pairs = set()
    if len(hostWords) > 0:
        if not isinstance(hostWords[-1],list):
            hostWords[-1] = [ hostWords[-1] ]
        tail = hostWords[-1]
        head = hostWords[:-1]
        for wordTail in tail:
            hw = head + [ wordTail ]
            hostNames = [ " ".join(hw[:endW]) for endW in range(len(hw),0,-1) ]
            for pair in zip(range(len(hostNames)),hostNames):
                pairs.add(pair)
    pairs = sorted(pairs)
    #print 'name: ', name, 'hostWords: ', hostWords, 'pairs: ', pairs
    return pairs

def splitNamesToHosts(names):
    """Take several full names, and return a single list of parts sorted by number of removed parts.
    The idea is to have more complete names to come first regardless of what full name they
    originated from.
    @param names list of full names, e.g. ['blackberry cv. apache','blackberry cv. apache str.234']
    @ret list of unique (num_of_removed_parts,prefix) pairs, e.g. [(0,'blackberry cv. apache'), (0,'blackberry cv. apache str.234')(1,'blackberry cv.'), (2, 'blackberry')]"""
    candNames = []
    for name in names:
        candNames.extend(splitNameToHosts(name))
    return sorted(set(candNames))

class VirHostParser:
    """Class that find phage-host pairs based on GenBank annotation as well as from phage and microbe names alone.
    The main entry point is assignHostsAndSave()"""

    def __init__(self,taxaTree):
        self.taxaTree = taxaTree
        taxaTree.makeNameIndex()
        self.virNode = taxaTree.getNode(virTaxid)
        self.topNode = taxaTree.getNode(virTaxid)
        self.unclassRootNode = taxaTree.getNode(unclassRootTaxid)
        #self.phRoot = taxaTree.getNode(phageTailedTaxid)

    def _excludeHost(self,node,hostNode):
        """Return True if host must be rejected.
        @param node viral node for which we are trying to find a host
        @param hostNode possible host node
        @todo this currently rejects giant mimivirus that is infected by a satellite virus"""
        return hostNode.isSubnode(self.virNode) or \
                hostNode.name.startswith('environmental') or \
                hostNode.isSubnode(self.unclassRootNode)

    def assignHost(self,node,candNames):
        """Search candidate names in the taxonomy tree.
        @param candNames [(edit_distance,name),(edit_distance,name),...]
        sorted by edit_distance.
        @ret list of nodes found.
        Going from smaller to larger edit distances, find all matches
        at a given edit distance and return nodes that are not descendants
        of any other found nodes. Once a list of nodes at a given edit distance
        is found, the method returns w/o looking at larger edit distances."""
        taxaTree = self.taxaTree
        groups = {}
        for pair in candNames:
            try:
                groups[pair[0]].append(pair[1])
            except KeyError:
                groups[pair[0]]=[pair[1]]
        searchMethods = (taxaTree.searchName,taxaTree.searchNames)
        for score in sorted(groups.keys()):
            names = groups[score]
            nodes = []
            for name in names:
                for method in searchMethods:
                    foundNodes = [ (score,foundNode) \
                        for foundNode in method(name) \
                        if not self._excludeHost(node,foundNode) ]
                    if len(foundNodes) > 0:
                        nodes.extend(foundNodes)
                        break
            if len(nodes) > 0:
                #@todo we should probably select the lower node if it is the 
                #only descendant (excluding the "unclassified" one) of the parent node.
                ns = [ node[1] for node in nodes ]
                ns = set(selectTopNodes(ns))
                nodes = [ node for node in nodes if node[1] in ns ]
                nodes = unique(nodes,stable=True)
                return nodes
        return []


    def assignHostOld(self,node,candNames):
        virNode = self.virNode
        virName = node.name
        #print "Node: %s" % virName
        for candName in candNames:
            #print "    %s : %s" % (virName, candName)
            miNodes = [ (candName[0],foundNode) for foundNode in self.taxaTree.searchName(candName[1]) \
                    if not self._excludeHost(node,foundNode)]
            if len(miNodes) > 0:
                #for miNode in miNodes:
                #    print "        miNode: ", miNode.name
                #break
                return miNodes
            miNodes = [ (candName[0],foundNode) for foundNode in self.taxaTree.searchNames(candName[1]) \
                    if not self._excludeHost(node,foundNode)]
            if len(miNodes) > 0:
                #for miNode in miNodes:
                #    print "        miNode: ", miNode.names
                #break
                return miNodes
        return []
    
    def assignHostsByVirName(self):
        for node in self.taxaTree.iterDepthTop(virTaxid):
            if not node.name.lower().startswith('environmental'):
                hostNodes = self.assignHost(node,candNames=splitNameToHosts(node.name))
                #print [ n[1].name for n in hostNodes ]
                if len(hostNodes) > 0:
                    if not hasattr(node,'host'):
                        node.host = Struct(nm=[])
                    node.host.nm = hostNodes

    def assignHostsByGbFeat(self,gbFile):
        taxaTree = self.taxaTree
        inGb = openCompressed(gbFile,'r')
        iRec = 0
        iHostRec = 0
        iFoundRec = 0
        for rec in SeqIO.parse(inGb,"genbank"):
            #if iFoundRec > 20:
            #    break
            #print "Rec start"
            candNames = []
            taxid = None
            for feat in rec.features:
                if feat.type == "source":
                    #print "Feat start"
                    #print str(feat)
                    #print feat.qualifiers.items()
                    quals = feat.qualifiers
                    if quals.has_key('specific_host'):
                        candNames.extend(quals['specific_host'])
                    if quals.has_key('lab_host'):
                        candNames.extend(quals['lab_host'])
                    if quals.has_key('host'):
                        candNames.extend(quals['host'])
                    if quals.has_key('db_xref'):
                        db_xrefs = quals['db_xref']
                        for db_xref in db_xrefs:
                            parts = db_xref.split('taxon:')
                            if len(parts) > 1:
                                taxid = int(parts[1])
            candNames = unique(candNames,stable=True)
            if len(candNames) > 0:
                iHostRec += 1
                #gi = int(rec.annotations['gi'])
                node = taxaTree.getNode(taxid)
                candNamesSplit = splitNamesToHosts(candNames)
                hostNodes = self.assignHost(node,candNamesSplit)
                if len(hostNodes) <= 0:
                    print "No matching host node for: taxid:%s hosts:%s" % (taxid,candNamesSplit)
                else:
                    if not hasattr(node,'host'):
                        node.host = Struct(gb=[])
                    if not hasattr(node.host,'gb'):
                        node.host.gb = []
                    #it is possible to have several GB records for one taxa
                    #we should actually handle situations like (1,node),(0,parent)
                    #coming from different GB records, but currently we just assume
                    #(n,node),(n,node)
                    node.host.gb = sorted(set(node.host.gb) | set(hostNodes))
                    iFoundRec += 1
            iRec += 1
            if iRec % 1 == 0:
                print "Processed %s GenBank records, found %s records with hosts, found host nodes for %s records" % (iRec,iHostRec,iFoundRec)
        inGb.close()

    def assignHosts(self,gbFile):
        self.assignHostsByGbFeat(gbFile=gbFile)
        self.assignHostsByVirName()

    def assignHostsAndSave(self,gbFile,outFile):
        self.assignHosts(gbFile=gbFile)
        hosts = {}
        for node in self.taxaTree.iterDepthTop(virTaxid):
            if hasattr(node,'host'):
                host = node.host
                if hasattr(host,'gb'):
                    hosts[node.id] = [ Struct(src='annot',taxid=hrec[1].id,score=hrec[0]) for hrec in host.gb ]
                elif hasattr(host,'nm'):
                    hosts[node.id] = [ Struct(src='name',taxid=hrec[1].id,score=hrec[0]) for hrec in host.nm ]
        dumpObj(hosts,outFile)


def loadHosts(inpFile,taxaTree):
    hosts = loadObj(inpFile)
    for idnode,hrecs in hosts.iteritems():
        try:
            node = taxaTree.getNode(idnode)
        except KeyError:
            print "Missing taxid in TaxaTree: %s" % idnode
        else:
            node.hosts = [ Struct(src=hrec.src,node=taxaTree.getNode(hrec.taxid),score=hrec.score) \
                    for hrec in hrecs ]
    return hosts 


# Methods related to selection of balanced sets of host-phage pairs
# and sequence sub-sampling


def _pickHostPhageSeqsOld(nodes,maxSeq=1):
    seqs = set()
    for node in nodes:
        # re-use already picked seqs wherever possible
        nodeSeqs = set([ (node.id,s.gi) for s in node.seq ])
        nodePickedSeqs = seqs & nodeSeqs
        nToPick = maxSeq - len(nodePickedSeqs)
        if nToPick > 0:
            nodeUnpickedSeqs = nodeSeqs-nodePickedSeqs
            seqSel = random.sample(nodeUnpickedSeqs,min(nToPick,len(nodeUnpickedSeqs)))
            seqs |= set(seqSel)
            assert len(seqSel) > 0
        else:
            print "DEBUG: All picked already"
    return seqs

def _pickHostPhageSeqs(node,maxSeq=1):
    if not hasattr(node,"seqPick"):
        node.seqPick = random.sample([ s.gi for s in node.seq ],min(maxSeq,len(node.seq)))

class PhageHostSeqPicker:
    """From all known known phage-host associations, pick a a stratified set of pairs.
    Currently, we pick one phage-host pair per microbial genus, wherever available.
    This is done to have a reasonable hope that a given phage has only a single true microbial host in the resulting set.
    Thus, all assignments by whatever classification algorithm to a different host can be considered as false positives.
    This class has two modes of operation. First one selects host/virus sequence pairs and can save them.
    Second one can load saved pairs into memory."""

    def __init__(self,taxaTree):
        """Constructor.
        @param taxaTree - taxonomy tree
        @pre sequences must be assigned by mapFastaRecordsToTaxaTree (unless load() will be called)
        @pre host data must be assigned by loadHosts (unless load() will be called)
        """
        self.taxaTree = taxaTree


    def pickSeqHosts(self):
        taxaTree = self.taxaTree
        micNodes = taxaTree.getNodes(micTaxids)
        micHosts = {} #dict with a set of phages for each host
        for node in taxaTree.iterDepthTop(virTaxid):
            #only viral nodes with known hosts and with RefSeq
            if hasattr(node,'hosts') and hasattr(node,'seq'):
                for hostRec in node.hosts:
                    #only annotated host assignment (not from name alone)
                    #and only microbial hosts
                    if hostRec.src == 'annot' and \
                            hostRec.node.isSubnodeAny(micNodes):
                        try:
                            micHosts[hostRec.node].add(node)
                        except KeyError:
                            micHosts[hostRec.node] = set([node])

        for host in micHosts:
            hostSpe = host.findRankInLineage("species")
            if hostSpe is not None:
                hostGen = host.findRankInLineage("genus")
                if hostGen is not None:
                    hostTop = hostGen
                else:
                    hostTop = hostSpe
                hostSeq = []
                # the node iterator will look at the host node itself first,
                # the in its subnodes, and then climb up to the genus subtree
                for hostNb in host.iterFillUp(topNode=hostTop):
                    if hasattr(hostNb,'seq'):
                        if hostNb is host or hostNb.isSubnode(host):
                            lowSeq = True
                        else:
                            lowSeq = False
                        #construct list of tuples such that reverse sorting it later
                        #will put the desired records (if any) first: those below
                        #the given host node and previously selected among them
                        hostSeq.append((lowSeq,hasattr(hostNb,'isSeqSel'),hostNb))
                hostSeq.sort(reverse=True)
                if len(hostSeq) > 0:
                    hostSeqSel = hostSeq[0][-1]
                    hostSeqSel.isSeqSel = True
                    host.nodeSeq = hostSeqSel
                    #host.nodeSpe = hostSpe
                    #host.nodeGen = hostGen

        seqHosts = {} #dict with a set of viruses for each selected host sequence
        for (host,vir) in micHosts.iteritems():
            if hasattr(host,'nodeSeq'):
                try:
                    seqHosts[host.nodeSeq] |= vir
                except KeyError:
                    seqHosts[host.nodeSeq] = vir.copy()
        self.seqHosts = seqHosts

    def groupSeqHosts(self):
        """Group seqHosts by genera, each value contains dict {host:set(vir)}"""
        seqHosts = self.seqHosts
        groups = {} 
        groupRanks = ("genus","subgenus","species")
        for host in seqHosts:
            groupNode = host.findLeftRankInLineage(groupRanks)
            assert groupNode is not None
            try:
                groups[groupNode][host] = seqHosts[host]
            except KeyError:
                groups[groupNode]= {host:seqHosts[host]}
        self.groups = groups

    def seqHostVirPairs(self):
        """Return a set([(vir,seqhost),...]) for a quick lookup of a given vir,host pair"""
        seqHosts = self.seqHosts
        pairs = set( ( (host,vir) for host in seqHosts for vir in seqHosts[host] ) )
        return pairs

    def seqHostPairsToMap(self,pairs):
        map = {}
        for pair in pairs:
            try:
                map[pair[0]].add(pair[1])
            except KeyError:
                map[pair[0]] = set([ pair[1] ])
        return map

    def seqHostVirs(self):
        return self.seqHosts

    def seqVirHosts(self):
        x = {}
        for (host,virs) in self.seqHosts.iteritems():
            for vir in virs:
                try:
                    x[vir].append(host)
                except KeyError:
                    x[vir] = [ host ]
        return x

    def seqVirHostPicks(self):
        return groupPairs(self.picks,keyField=1)

    def seqHostVirPicks(self):
        return groupPairs(self.picks,keyField=0)

    def printGroupSeqHosts(self):
        """Print dict of dicts of lists {mic genus : mic species : phages}"""
        groups = self.groups
        for (miGen,miSpes) in groups.items():
            print "%s %s" % (miGen.name, miGen.rank)
            for (miSpe,virs) in miSpes.items():
                print "\t%s %s" % (miSpe.name,miSpe.rank)
                for vir in virs:
                    print "\t\t%s" % vir.name

    def pickPairs(self,maxMicSpe=1,maxMicSeq=1,maxVir=1,maxVirSeq=1):
        groups = self.groups
        #mics = set()
        #virs = set()
        picks = set()
        for (miGen,miSpes) in groups.items():
            print "%s %s" % (miGen.name, miGen.rank)
            miSpeNodes = miSpes.keys()
            miSpeSel = random.sample(miSpeNodes,min(maxMicSpe,len(miSpeNodes)))
            #mics |= set(miSpeSel)
            for miSpe in miSpeSel:
                _pickHostPhageSeqs(node=miSpe,maxSeq=maxMicSeq)
                viNodes = miSpes[miSpe]
                viSel = random.sample(viNodes,min(maxVir,len(viNodes)))
                #virs |= set(viSel)
                print "\t%s %s" % (miSpe.name,miSpe.rank)
                for vir in viSel:
                    print "\t\t%s" % vir.name
                    _pickHostPhageSeqs(node=vir,maxSeq=maxVirSeq)
                    picks.add((miSpe,vir))
        self.picks = picks

    def checkPairs(self,giToTaxa):
        taxaTree = self.taxaTree
        picks = self.picks
        mics = set( ( x[0] for x in picks ) )
        virs = set( ( x[1] for x in picks ) )
        for pickSet in (mics,virs):
            for node in pickSet:
                for gi in node.seqPick:
                    assert node.id == giToTaxa[gi]
        assert mics.issubset(set(self.seqHostVirs())) and virs.issubset(set(self.seqVirHosts()))

    def save(self,out):
        picks = self.picks
        picks = [ ((mic.id,mic.seqPick),(vir.id,vir.seqPick)) for (mic,vir) in picks ]
        pairs = self.seqHostVirPairs()
        pairs = [ (mic.id,vir.id) for (mic,vir) in pairs ]
        dumpObj(Struct(picks=picks,pairs=pairs),out)

    def load(self,inp):
        data = loadObj(inp)
        taxaTree = self.taxaTree
        pairs = data.pairs
        pairs = [ (taxaTree.getNode(micid), taxaTree.getNode(virid)) for (micid,virid) in pairs ]
        picks = []
        for pick in data.picks:
            mic = taxaTree.getNode(pick[0][0])
            mic.seqPick = pick[0][1]
            vir = taxaTree.getNode(pick[1][0])
            vir.seqPick = pick[1][1]
            picks.append((mic,vir))
        self.picks = picks
        self.seqHosts = self.seqHostPairsToMap(pairs=pairs)

    def saveSeqIds(self,out):
        """Save a set of all (microbial and viral) sequence IDs (GIs) selected by pickPairs().
        This is for pulling all relevant sequences into internal representation."""
        picks = self.picks
        mics = set( ( x[0] for x in picks ) )
        virs = set( ( x[1] for x in picks ) )
        nodes = mics | virs
        seqIds = []
        for node in nodes:
            seqIds += node.seqPick
        dumpObj(set(seqIds),out)

    def makeIdLabsPicks(self,div="all"):
        """Return IdLables(gi->taxid) object for microbial, viral or all sequences selected by pickPairs()"""
        picks = self.picks
        mics = set( ( x[0] for x in picks ) )
        virs = set( ( x[1] for x in picks ) )
        if div == "mic":
            nodes = mics
        elif div == "vir":
            nodes = virs
        elif div == "all":
            nodes = mics | virs
        else:
            raise ValueError(div)
        idToTaxid = []
        for node in nodes:
            idToTaxid += [ (idSeq,node.id,0) for idSeq in node.seqPick ]
        return IdLabels(records=idToTaxid)


