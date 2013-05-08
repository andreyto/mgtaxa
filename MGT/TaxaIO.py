### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.TaxaTree import TaxaNode
import json

class NodeStorageNcbiDump:
    """Loads tree nodes from NCBI flat dump file(s).
    This is 3x slower than loading through NodeStoragePickle representation."""

    fields = \
    (
        ('taxid',int),
        ('partaxid',int),
        ('rank',str),
        ('embl_code',str),
        ('divid',int),
#        ('inh_div',bool),
#        ('gcode_id',int),
#        ('inh_gc',bool),
#        ('mgcode_id',int),
#        ('inhmgc',bool),
#        ('gbhidden',bool),
#        ('hidsubtree',bool),
#        ('comments',str)
     )

    nInputFields = 5

    # We will be parsing this string:
    # '1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n'
    
    def __init__(self,ncbiDumpFile=None,ncbiNamesDumpFile=None,ncbiMergedDumpFile=None,allNames=False):
        """Load taxonomy tree nodes from NCBI dump files.
        @param ncbiDumpFile nodes.dump
        @param ncbiNamesDumpFile names.dump
        @param allNames whether to load all names to node attribute 'names'"""
        self.ncbiDumpFile = ncbiDumpFile
        self.ncbiNamesDumpFile = ncbiNamesDumpFile
        self.ncbiMergedDumpFile = ncbiMergedDumpFile
        self.allNames = allNames
    
    def _loadNodes(self,ncbiDumpFile):
        inp = openCompressed(ncbiDumpFile,'r')
        nodes = {}
        n_splits = self.nInputFields
        #delimRe = re.compile(r"\s*\|\s*")
        for rec in inp:
            # On profiling, string methods outperformed regexes (compare next two lines):
            #values = [ x.strip() for x in rec.split('\t|\t',n_splits)[:n_splits]]
            #values = delimRe.split(rec,n_splits)
            
            # This will not split the last field, but we never use it anyway:
            values = rec.split("\t|\t")
            node = TaxaNode()
            node.id = int(values[0])
            node.idpar = int(values[1])
            node.rank = values[2].replace(' ','_')
            node.divid = int(values[4])
            #in NCBI file, root node points to itself as a parent.
            #We replace it with 0 for consistency with our SQL DB representation, where
            #circular self-reference would be inconvenient.
            if node.idpar == node.id:
                node.idpar = 0
            nodes[node.id] = node
            assert node.id < ncbiTaxidMax, "We assume that dump file is pristine NCBI file and "+\
                    "assert that taxonomy ID (%s) < our max limit (%s)" % (node.id,ncbiTaxidMax)
        inp.close()
        self.nodes = nodes

    def _loadNames(self,ncbiNamesDumpFile):
        nodes = self.nodes
        inp = openCompressed(ncbiNamesDumpFile,'r')
        allNames = self.allNames
        for line in inp:
            #rec = [ x.strip() for x in line.split('|') ]
            rec = line.split("\t|\t")
            if rec[3].startswith("scientific name"):
                nodes[int(rec[0])].name = rec[1]
            if allNames:
                node = nodes[int(rec[0])]
                try:
                    node.names[rec[1]] = rec[3]
                except AttributeError:
                    node.names = { rec[1] : rec[3] }
        inp.close()

    def _loadMerged(self,ncbiMergedDumpFile):
        inp = openCompressed(ncbiMergedDumpFile,'r')
        data = dict(( (int(rec[0].strip()),int(rec[1].strip())) for rec in \
                ( line.split("|") for line in inp ) ))
        inp.close()
        return data

    def load(self):
        self._loadNodes(self.ncbiDumpFile)
        if self.ncbiNamesDumpFile is not None:
            self._loadNames(self.ncbiNamesDumpFile)
        if self.ncbiMergedDumpFile is not None:
            merged = self._loadMerged(self.ncbiMergedDumpFile)
        else:
            merged = {}
        return dict(nodes=self.nodes,merged=merged)


class NodeStorageNewick:
    import string

    def __init__(self,fileName,labeler=lambda node: node.name):
        self.fileName = fileName
        self.labeler = labeler
        trans = list(allChr())
        allowedChr = set(string.letters + string.digits)
        for i in range(len(trans)):
            if trans[i] not in allowedChr:
                trans[i] = "_"
        self.trans = ''.join(trans)

    def _quote(self,label):
        return "'%s'" % label.translate(self.trans).replace("'","''")
        
    def save(self,tree):
        labeler = self.labeler
        out = openCompressed(self.fileName,"w")
        lastVisit = 1
        for (node,visit) in tree.iterDepthTwice():
            if visit == 1:
                if lastVisit == 2:
                    out.write(",")
            if node.isLeaf():
                if visit == 2:
                    out.write(self._quote(labeler(node)))
            else:
                if visit == 1:
                    out.write("(")
                else:
                    out.write(")%s\n" % self._quote(labeler(node)))
            lastVisit = visit
        out.close()


class NodeStorageHypView:
    import string

    def __init__(self,fileName,labeler=lambda node: node.name):
        self.fileName = fileName
        self.labeler = labeler
        trans = list(allChr())
        allowedChr = set(string.letters + string.digits)
        for i in range(len(trans)):
            if trans[i] not in allowedChr:
                trans[i] = "_"
        self.trans = ''.join(trans)

    def _quote(self,label):
        return "'%s'" % label.translate(self.trans).replace("'","''")
        
    def save(self,tree):
        labeler = self.labeler
        out = openCompressed(self.fileName,"w")
        for node in tree.iterDepthTop():
            out.write("%s %s 0 html\n" % (node.getDepth(),self._quote(labeler(node))))
        out.close()
    

class NodeStoragePickleDict:
    """
    Implements NodeStorage interface through Python pickle mechanism.
    This pickles entire node dict"""

    def __init__(self,fileName):
        self.fileName = fileName

    def save(self,tree):
        nodes = tree.getNodesDict()
        merged = tree.getMerged()
        # TaxaTreeNode.__getstate__() now takes care of skipping node cross-references
        #for node in nodes.itervalues():
        #    del (node.par, node.children)
        dumpObj(dict(nodes=nodes,merged=merged),self.fileName)

    def load(self):
        data = loadObj(self.fileName)
        nodes = data["nodes"]
        for node in nodes.itervalues():
            node.children = []
        return data


class NodeStoragePickle:
    """
    Implements NodeStorage interface through Python pickle mechanism.
    This pickles nodes one-by-one"""

    def __init__(self,fileName):
        self.fileName = fileName

    def save(self,tree):
        nodes = tree.getNodesDict()
        merged = tree.getMerged()
        # TaxaTreeNode.__getstate__() now takes care of skipping node cross-references
        #for node in nodes.itervalues():
        #    del (node.par, node.children)
        out = openCompressed(self.fileName,"w")
        dump(merged,out,-1)
        for node in nodes.itervalues():
            dump(node,out,-1)
        out.close()

    def load(self):
        #nodes = loadObj(self.fileName)
        nodes = {}
        inp = openCompressed(self.fileName,"rb")
        merged = load(inp)
        try:
            while True:
                node = load(inp)
                node.children = []
                nodes[node.id] = node
        except:
            pass
        inp.close()
        return dict(nodes=nodes,merged=merged)

##@todo there are tons of Python wrappers around streaming (SAX-like) JSON C libs,
##e.g. http://pypi.python.org/pypi/ijson/. They are supposed to be much faster,
##and the underlying C lib is better supported than jsonlib

class NodeStorageJsonPy:
    """
    Implements NodeStorage interface through json module serialization.
    This pickles nodes as list"""

    def __init__(self,fileName):
        self.fileName = fileName

    def save(self,tree):
        import json
        nodes = tree.getNodesDict()
        merged = tree.getMerged()
        # TaxaTreeNode.__getstate__() now takes care of skipping node cross-references
        #for node in nodes.itervalues():
        #    del (node.par, node.children)
        out = openCompressed(self.fileName,"w")
        #json saves keys as strings, so there would be little point in
        #dumping the nodes dict or merged dict
        json.dump(dict(nodes=nodes.values(),merged=merged.items()),
                out,
                check_circular=False,
                separators=(',', ':'),
                default=TaxaNode.__getstate__)
        out.close()

    def load(self):
        import json
        inp = openCompressed(self.fileName,"rb")
        data = json.load(inp,object_hook=lambda o: TaxaNode(**o) if "idpar" in o else o)
        print "JsonDict Loaded the nodes"
        inp.close()
        data["nodes"] = dict( ((node.id,node) for node in data["nodes"]) )
        data["merged"] = dict( (item for item in data["merged"]) )
        print "Made dict"
        return data

class NodeStorageJsonUjson:
    """
    Implements NodeStorage interface through json module serialization.
    This pickles nodes as list"""

    def __init__(self,fileName):
        self.fileName = fileName

    def save(self,tree):
        import json
        nodes = tree.getNodesDict()
        merged = tree.getMerged()
        # TaxaTreeNode.__getstate__() now takes care of skipping node cross-references
        #for node in nodes.itervalues():
        #    del (node.par, node.children)
        out = openCompressed(self.fileName,"w")
        #json saves keys as strings, so there would be little point in
        #dumping the nodes dict
        json.dump(dict(nodes=nodes.values(),merged=merged.items()),
                out,check_circular=False,
                separators=(',', ':'),
                default=TaxaNode.__getstate__)
        out.close()

    def load(self):
        import ujson
        inp = openCompressed(self.fileName,"rb")
        buf = inp.read()
        inp.close()
        #nodes = json.loads(buf,object_hook=lambda o: TaxaNode(**o) if "idpar" in o else o)
        #return dict( ((node.id,node) for node in nodes) )
        data = ujson.loads(buf)
        data["nodes"] = dict( ((node["id"],TaxaNode(**node)) for node in data["nodes"]) )
        data["merged"] = dict( (item for item in data["merged"]) )
        return data

class NodeStorageJsonLib:
    """
    Implements NodeStorage interface through json module serialization.
    This pickles nodes as list"""

    def __init__(self,fileName):
        self.fileName = fileName

    def save(self,tree):
        import json
        nodes = tree.getNodesDict()
        merged = tree.getMerged()
        # TaxaTreeNode.__getstate__() now takes care of skipping node cross-references
        #for node in nodes.itervalues():
        #    del (node.par, node.children)
        out = openCompressed(self.fileName,"w")
        #json saves keys as strings, so there would be little point in
        #dumping the nodes dict
        json.dump(dict(nodes=nodes.values(),merged=merged.items()),
                out,check_circular=False,
                separators=(',', ':'),
                default=TaxaNode.__getstate__)
        out.close()

    def load(self):
        import jsonlib
        inp = openCompressed(self.fileName,"rb")
        buf = inp.read()
        inp.close()
        #nodes = json.loads(buf,object_hook=lambda o: TaxaNode(**o) if "idpar" in o else o)
        #return dict( ((node.id,node) for node in nodes) )
        data = jsonlib.read(buf,use_float=True)
        data["nodes"] = dict( ((node["id"],TaxaNode(**node)) for node in data["nodes"]) )
        data["merged"] = dict( (item for item in data["merged"]) )
        return data

class NodeStorageJsonLines:
    """
    Implements NodeStorage interface through json module serialization.
    This pickles nodes one-by-one"""

    def __init__(self,fileName):
        self.fileName = fileName

    def save(self,tree):
        nodes = tree.getNodesDict()
        # TaxaTreeNode.__getstate__() now takes care of skipping node cross-references
        #for node in nodes.itervalues():
        #    del (node.par, node.children)
        out = openCompressed(self.fileName,"w")
        #json saves keys as strings, so there would be little point in
        #dumping the nodes dict
        for node in nodes.itervalues():
            json.dump(node.__getstate__(),out,check_circular=False,separators=(',', ':'))
            out.write("\n")
        out.close()

    def load(self):
        import jsonlib
        inp = openCompressed(self.fileName,"rb")
        #nodes = json.loads(buf,object_hook=lambda o: TaxaNode(**o) if "idpar" in o else o)
        #return dict( ((node.id,node) for node in nodes) )
        nodes = {}
        for line in inp:
            node = jsonlib.read(line.strip(),use_float=True)
            nodes[node["id"]] = TaxaNode(**node)
        inp.close()
        return nodes

NodeStorageJson = NodeStorageJsonPy

class LibSeaElement:

    def __init__(self,out):
        self.out = out

    def close(self):
        self.out.write(\
        dedent(\
        """
        ];
        """)
        )


class LibSeaGraph(LibSeaElement):

    def open(self,**kw):
        self.out.write(\
        dedent(\
        """
        Graph
        {
            ### metadata ###
            @name="%(name)s";
            @description="%(description)s";
            @numNodes=%(numNodes)i;
            @numLinks=%(numLinks)i;
            @numPaths=0;
            @numPathLinks=0;
            
            ### structural data ###
        """ % kw)
        )

    def close(self):
        self.out.write(\
        dedent(\
        """
        ### visualization hints ###
        @filters=;
        @selectors=;
        @displays=;
        @presentations=;
        
        ### interface hints ###
        @presentationMenus=;
        @displayMenus=;
        @selectorMenus=;
        @filterMenus=;
        @attributeMenus=;
        }
        """)
        )

class LibSeaLinks(LibSeaElement):

    def open(self):
        self.out.write(\
        dedent(\
        """
        @links=
        [
        """)
        )

    def link(self,comma,**kw):
        self.out.write(comma+"\n{ @source=%(source)i; @destination=%(destination)i; }" % kw)


class LibSeaPaths(LibSeaElement):
    
    def open(self):
        self.out.write("@paths=;")

    def close(self):
        pass
     

class LibSeaAttributes(LibSeaElement):

    def open(self):
        self.out.write(\
        dedent(\
        """
        ### attribute data ###
        @enumerations=;
        @attributeDefinitions=[
        """)
        )

class LibSeaRootAttribute(LibSeaElement):

    def open(self,**kw):
        self.out.write(\
        dedent(\
        """
        {
         @name=$root;
         @type=bool;
         @default=|| false ||;
         @nodeValues=[ { @id=%(id)i; @value=T; } ];
         @linkValues=;
         @pathValues=;
        }""" % kw)
        )

    def close(self):
        pass

class LibSeaTreeLinkAttribute(LibSeaElement):

    def open(self):
        self.out.write(\
        dedent(\
        """
        {
         @name=$tree_link;
         @type=bool;
         @default=|| false ||;
         @nodeValues=;
         @linkValues= [
        """)
        )

    def linkValue(self,comma,**kw):
        self.out.write(comma+"\n{ @id=%(id)i; @value=T; }" % kw)

    def close(self):
        self.out.write(\
        dedent(\
        """
         ];
         @pathValues=;
        }""")
        )


class LibSeaNodeNamesAttribute(LibSeaElement):

    def open(self):
        self.out.write(\
        dedent(\
        """
        {
         @name=$name_text;
         @type=string;
         @default=|| "unknown" ||;
         @nodeValues=[
        """)
        )

    def nodeValue(self,comma,**kw):
        self.out.write(comma+"""\n{ @id=%(id)i; @value="name=%(name)s"; }""" % kw)

    def close(self):
        self.out.write(\
        dedent(\
        """
         ];
         @linkValues=;
         @pathValues=;
        }""")
        )

class LibSeaQualifiers(LibSeaElement):

    def open(self):
        self.out.write(\
        dedent(\
        """
        @qualifiers=[
            {
                @type=$spanning_tree;
                @name=$sample_spanning_tree;
                @description=;
                @attributes=[
                    { @attribute=0; @alias=$root; },
                    { @attribute=1; @alias=$tree_link; }
                ];
            }
        ];
        """)
        )

    def close(self):
        pass


class NodeStorageLibSea:
    import string

    def __init__(self,fileName,labeler=lambda node: node.name):
        self.fileName = fileName
        self.labeler = labeler
        trans = list(allChr())
        allowedChr = set(string.letters + string.digits)
        for i in range(len(trans)):
            if trans[i] not in allowedChr:
                trans[i] = "_"
        self.trans = ''.join(trans)

    def _quote(self,label):
        return "'%s'" % label.translate(self.trans).replace("'","''")
        
    def save(self,tree):
        labeler = self.labeler
        out = openCompressed(self.fileName,"w")
        graphProp = {}
        graphProp["numNodes"] = tree.numNodes()
        graphProp["numLinks"] = tree.numLinks()
        graphProp["name"] = "Graph"
        graphProp["description"] = "Taxonomy"
        graph = LibSeaGraph(out)
        graph.open(**graphProp)
        links = LibSeaLinks(out)
        links.open()
        libSea_nodeId = 0
        comma = ""
        it = tree.iterDepthTop()
        node = it.next()
        node.libSea_nodeId = libSea_nodeId
        libSea_nodeId += 1
        for node in it:
            node.libSea_nodeId = libSea_nodeId
            links.link(comma,**{"source":node.getParent().libSea_nodeId, "destination":libSea_nodeId})
            comma = ","
            libSea_nodeId += 1
        links.close()
        paths = LibSeaPaths(out)
        paths.open()
        paths.close()
        attributes = LibSeaAttributes(out)
        attributes.open()
        rootAttribute = LibSeaRootAttribute(out)
        rootAttribute.open(id=0)
        rootAttribute.close()
        out.write(",\n")
        treeLinkAttribute = LibSeaTreeLinkAttribute(out)
        treeLinkAttribute.open()
        comma = ""
        for id in xrange(graphProp["numLinks"]):
            treeLinkAttribute.linkValue(comma,id=id)
            comma = ","
        treeLinkAttribute.close()
        out.write(",\n")
        nodeNamesAttribute = LibSeaNodeNamesAttribute(out)
        nodeNamesAttribute.open()
        comma = ""
        for node in tree.iterDepthTop():
            nodeNamesAttribute.nodeValue(comma,id=node.libSea_nodeId,name=self._quote(labeler(node)))
            comma = ","
        nodeNamesAttribute.close()
        attributes.close()
        qualifiers = LibSeaQualifiers(out)
        qualifiers.open()
        qualifiers.close()
        graph.close()
        out.close()
