from MGT.Common import *
from MGT.TaxaTree import TaxaNode

class NodeStorageNcbiDump:
    """@todo All attempts to create a faster loading representation for nodes were unsuccessful.
    cPickle, Numpy arrays as intermediate buffer for tree nodes were slower than parsing flat file."""


    fields = \
    (
        ('taxid',int),
        ('partaxid',int),
        ('rank',str),
#        ('embl_code',str),
#        ('divid',int),
#        ('inh_div',bool),
#        ('gcode_id',int),
#        ('inh_gc',bool),
#        ('mgcode_id',int),
#        ('inhmgc',bool),
#        ('gbhidden',bool),
#        ('hidsubtree',bool),
#        ('comments',str)
     )

    nInputFields = 3

    ## We will be parsing this string:
    ## '1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n'
    
    def __init__(self,ncbiDumpFile=None,ncbiNamesDumpFile=None):
        """Load taxonomy tree nodes from NCBI dump files"""
        self.ncbiDumpFile = ncbiDumpFile
        self.ncbiNamesDumpFile = ncbiNamesDumpFile
    
    def _loadNodes(self,ncbiDumpFile):
        inp = openCompressed(ncbiDumpFile,'r')
        nodes = {}
        n_splits = self.nInputFields
        #delimRe = re.compile(r"\s*\|\s*")
        for rec in inp:
            ## On profiling, string methods outperformed regexes (compare next two lines):
            #values = [ x.strip() for x in rec.split('\t|\t',n_splits)[:n_splits]]
            #values = delimRe.split(rec,n_splits)
            
            # This will not split the last field, but we never use it anyway:
            values = rec.split("\t|\t")
            node = TaxaNode()
            node.id = int(values[0])
            node.idpar = int(values[1])
            node.rank = values[2].replace(' ','_')
            #in NCBI file, root node points to itself as a parent.
            #We replace it with 0 for consistency with our SQL DB representation, where
            #circular self-reference would be inconvenient.
            if node.idpar == node.id:
                node.idpar = 0
            nodes[node.id] = node
        inp.close()
        self.nodes = nodes

    def _loadNames(self,ncbiNamesDumpFile):
        nodes = self.nodes
        inp = openCompressed(ncbiNamesDumpFile,'r')
        for line in inp:
            #rec = [ x.strip() for x in line.split('|') ]
            rec = line.split("\t|\t")
            if rec[3].startswith("scientific name"):
                nodes[int(rec[0])].name = rec[1]
        inp.close()

    def load(self):
        self._loadNodes(self.ncbiDumpFile)
        if self.ncbiNamesDumpFile is not None:
            self._loadNames(self.ncbiNamesDumpFile)
        return self.nodes


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
    

class LibSeeElement:

    def __init__(self,out):
        self.out = out

    def close(self):
        self.out.write(\
        dedent(\
        """
        ];
        """)
        )


class LibSeaGraph(LibSeeElement):

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

class LibSeaLinks(LibSeeElement):

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


class LibSeaPaths(LibSeeElement):
    
    def open(self):
        self.out.write("@paths=;")

    def close(self):
        pass
     

class LibSeaAttributes(LibSeeElement):

    def open(self):
        self.out.write(\
        dedent(\
        """
        ### attribute data ###
        @enumerations=;
        @attributeDefinitions=[
        """)
        )

class LibSeaRootAttribute(LibSeeElement):

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

class LibSeaTreeLinkAttribute(LibSeeElement):

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


class LibSeaNodeNamesAttribute(LibSeeElement):

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

class LibSeaQualifiers(LibSeeElement):

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
