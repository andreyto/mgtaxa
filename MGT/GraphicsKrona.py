"""Generator of input files for Krona JavaScript taxonomic charting library.
Krona citation:
Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385
Krona code:
http://krona.sourceforge.net
"""

from MGT.Common import *


class KronaWriter(object):
    """Class to write Krona XML input"""
    _kronaEnterTpl="""\
<krona collapse="true" key="false">
<attributes magnitude="count">
 <attribute display="Magnitude">magnitude</attribute>
 <attribute display="Count">count</attribute>
 <attribute display="Rank">rank</attribute>
 <attribute>score</attribute>
 <attribute display="Tax ID" hrefBase="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=" target="taxonomy">taxid</attribute>
</attributes>
<datasets>
  <dataset>MGTAXA taxonomic assignments</dataset>
</datasets>
"""

    #<attribute display="Avg. confidence">score</attribute>
    #<color attribute="score" hueStart="0" hueEnd="120" valueStart="0.413" valueEnd="0.932367571822115" default="false" ></color>

    _kronaExitTpl = """\
</krona>
"""

    _nodeEnterTpl="""\
    <node name="%(name)s">
    <count><val>%(count)s</val></count>
    <score><val>%(score)s</val></score>
    <rank><val>%(rank)s</val></rank>
    <taxid><val href="%(taxid_href)s">%(taxid)s</val></taxid>
    """
    _nodeExitTpl="""\
    </node>
    """

    def __init__(self,taxaTree):
        self.taxaTree = taxaTree

    def getTaxaTree(self):
        return self.taxaTree

    def initCount(self):
        taxaTree = self.taxaTree
        taxaTree.setAttribute("mgt_cnt",0)

    def addSamples(self,samples):
        for sample in samples:
            self.addSample(sample=sample)

    def addSample(self,sample):
            taxid = sample[1]
            weight = sample[2]
            if taxid >= validTaxid:
                node = self.taxaTree.getNode(taxid)
                node.mgt_cnt += weight

    def finishCount(self):
        taxaTree = self.getTaxaTree()
        root = taxaTree.getRootNode()
        root.setTotal("mgt_cnt","mgt_cnt_t")

    def initOut(self,xmlOut):
        self.xmlOut = xmlOut
        self.out = openCompressed(xmlOut,"w")
        self.out.write(self._kronaEnterTpl)

    def writeNodes(self):
        out = self.out
        taxaTree = self.getTaxaTree()
        root = taxaTree.getRootNode()
        
        def _nodeWriter(node,visit,context=self):
            #if no total count here, write nothing and signal not to go down as well
            if not node.mgt_cnt_t:
                return True
            #skip noRanks but go down (flatten, because they contributed to upper nodes already)
            if node.rank != noRank: 
                if visit==1:
                    context.out.write(context._nodeEnterTpl % dict(name=node.name,
                            count=node.mgt_cnt_t,
                            score=0,
                            rank=node.rank,
                            taxid_href=node.id if \
                                    node.id < ncbiTaxidMax and node.id >= validTaxid \
                                    else \
                                    rootTaxid,
                            taxid=node.id))
                else:
                    context.out.write(context._nodeExitTpl)
        
        root.visitDepthTwice(_nodeWriter)


    def finishOut(self):
        self.out.write(self._kronaExitTpl)
        self.out.close()

    def cleanTree(self):
        self.taxaTree.getRootNode().delAttributes(("mgt_cnt","mgt_cnt_t"))

    def genKrona(self,htmlOut):
        """Generate final Krona file"""
        cmd = options.krona.xmlToChart+["-u","/static/krona","-o",htmlOut,self.xmlOut]
        run(cmd,debug=True)

    def write(self,htmlOut):
        self.initOut(xmlOut=htmlOut+".krona.xml")
        self.writeNodes()
        self.finishOut()
        #self.cleanTree()
        self.genKrona(htmlOut=htmlOut)
