from MGT.JCVI.MgxAnnotFormat import *
from MGT.GFF import GFF3Record, GFF3Attributes, GFF3Header
from MGT.GFFTools import GFF3Graphics
from MGT.Taxa import *
from MGT.ImmClassifierApp import TaxaPred


def gffBtab(inpBtab,predTaxa,taxaTree,predMinLenSamp,outDir):

    inpBtab = BtabReader(inpBtab)
    idSampToPred = indexTuples(it.izip(predTaxa.idSamp,predTaxa.predTaxid,predTaxa.lenSamp))
    lastIdSamp = ""
    lastIdQ = ""
    outGff = None
    recGff = None
    for rec in inpBtab:
        id_q = rec.parsed["id_q"].strip()
        orf_cont = parseOrfOnContigId(id_q)
        idSamp = orf_cont["id_cont"]
        if lastIdQ == id_q:
            continue
        lastIdQ = id_q
        if idSamp in idSampToPred:
            pred = idSampToPred[idSamp]
            if pred[2] >= predMinLenSamp:
                if lastIdSamp != idSamp:
                    lastIdSamp = idSamp
                    if outGff:
                        outGff.close()
                    taxid = pred[1]
                    taxaName = taxaTree.getNode(taxid).name if taxaTree and taxid > 0 else "unassigned"
                    outFileRt = pjoin(outDir,"%s.%s" % \
                            idSamp,
                            taxaName.replace(" ","_").replace("/","_"))
                    gffFile = outFileRt + ".gff3"
                    if os.path.exists(gffFile):
                        outGff = open(gffFile,"a")
                    else:
                        outGff = open(gffFile,"w")
                        outGff.write(str(GFF3Header()))
                    recGff = GFF3Record(seqid=idSamp)
                recGff.type = "protein_match"
                recGff.start = orf_cont["start_orf"]
                recGff.end = orf_cont["end_orf"]
                recGff.strand = '+' if orf_cont["strand_orf"] > 0 \
                        else '-' if orf_cont["strand_orf"] < 0 \
                        else '.'
                # no attributes should be automatically carried forward from the previous values:
                recGff.attribs = GFF3Attributes()
                ats = recGff.attribs
                parsed = rec.parse()
                ats["ID"] = parsed["id_s"]
                ats["Name"] = parsed["descr"]
                if len(ats["Name"]) == 0:
                    ats["Name"] = ats["ID"]
                outGff.write(str(recGff))

    if outGff:
        outGff.close()


class AnnotTaxaClassifier(object):
    
    #weights per hit confidence type to affect voting
    confMapper = {
            "Reviewed":1.0,
            "HighConfidence":0.9,
            "Putative":0.8,
            "ConservedDomain":0.7,
            "LowConfidence":0.5,
            "LowestConfidence":0.25
            }

    def __init__(self,taxaTree,taxaLevels):
        self.taxaTree = taxaTree
        self.taxaLevels = taxaLevels
        self.nameStats = defdict(int)

    def getNode(self,name):
        """Get uniqie tree node from the taxa name"""
        annNodes = self.taxaTree.searchName(name,fuzzy=True)
        #print "DEBUG: tax_name=[%s] node_names=%s" % (name,
        #        [(node.rank,node.name) for node in annNodes])
        self.nameStats[len(annNodes)]+=1
        if len(annNodes) == 1:
            return annNodes[0]
        else:
            return None

    def predict(self,annotRecs):
        taxaLevels = self.taxaLevels
        linNodesConf = defdict(float)
        linNodesCnt = defdict(int)
        annotNodes = {}
        cntHits = 0
        confHits = 0
        for rec in annotRecs:
            node = self.getNode(rec["tax_name"])
            annotNodes[id(rec)] = node
            conf_x = self.confMapper[rec["conf"]]
            if node:
                skip = False
                #there are hits to environmental nodes directly under tree root,
                #those can dillute the taxonomically informative assignments and
                #should be skipped
                if not node.findRankInLineage("superkingdom"):
                    skip = True
                if not skip:
                    confHits += conf_x
                    cntHits += 1
                    for linNode in node.lineage():
                        linNodesConf[linNode] += conf_x
                        linNodesCnt[linNode] += 1
        linNodesConf = dict(linNodesConf)
        linNodesCnt = dict(linNodesCnt)
        #sort by (-rank,cnt)
        sortedLin = sorted([(taxaLevels.getLevelId(node.linn_level),-float(confHitsNode)/confHits,node) \
                for (node,confHitsNode) in linNodesConf.items()])
        ruleVars=dict(cntHitsNode=0,
            score=0,
            confHits=confHits,
            cntHits=cntHits)
        ret = dict(annotNodes=annotNodes,predNode=None,ruleVars=ruleVars)
        for (idLev,score,node) in sortedLin:
            assert idLev >= taxaLevels.minLinnId
            score = - score
            cntHitsNode = linNodesCnt[node]
            ruleVars["score"] = score
            ruleVars["cntHitsNode"] = cntHitsNode
            if cntHits == cntHitsNode or \
                    (cntHits >= 4 and cntHitsNode >= 2 and score >= 0.7499999999):
                ret["predNode"] = node
                break
        return ret

    def printStats(self):
        print "DEBUG: nameStats=%s" % (sorted(self.nameStats.items()),)

class AnnotTaxaFileClassifier(object):

    _annotTaxaFileClassifierDoc = """\
General description
-------------------
-------------------

This directory contains files that are created by processing the output from JCVI metagenomic 
annotation pipeline and, optionally, MGTAXA taxonomic predictions for long assembly contigs.
Analysis of the annotation output implements a majority voting scheme for BLASTP hits in 
order to assign a consensus taxonomy to a contig with multiple predicted ORFs.

Method
------
------

The voting scheme works like this, for a given contig:
1. Pick one best BLASTP hit for each ORF from camera_annotation_parser.raw.combined.out.gz 
pipeline output file, according to JCVI confidence assignments (such as HighConfidence, 
ConservedDomain etc designations - they already take into an account identity, coverage and score of a hit).
2. Find a node in NCBI taxonomic tree that has the lowest rank and that is in the lineage of 
at least 75% of ORF hits from p.1. A fixed value of 75% is something that merely seemed 
reasonable; there are could be other choices. When those 75% are computed, each ORF contributes 
with a weight adjusted according to its confidence value assigned by the pipeline. The weights 
are arbitrary parameters; those were used so far:
    	confMapper = {
            	"Reviewed":1.0,
            	"HighConfidence":0.9,
            	"Putative":0.8,
            	"ConservedDomain":0.7,
            	"LowConfidence":0.5,
            	"LowestConfidence":0.25
            	}

There is an extra requirement that there are at least four ORFs in the contig, and at least two of 
them contributed to the taxonomic node that is being considered. Otherwise, the contig is assigned
based on a total consensus of all ORFs (100% vote).

Input
-----
-----
It reads various files generated by the JCVI annotation pipeline, as well as the contig FASTA that 
served as input to the pipeline (to obtain the contig lengths). To avoid creating large hashes, it
has pre-conditions on the order of records in the annotation pipeline outputs. Specifically, it is
expected that the records in the annotation files are in the original order of the input contigs.
The program does keep a hash of contig names seen so far in the annotation files and should exit
with an error if a discontinuous set of records for a given contig is detected.

Output
------
------

The main output file is called, by default, annot.csv
The file contains predictions and the selected hit records that contributed to a decision, 
each with full main-rank lineage.

Because there are different types of records, the format is not a set of columns with fixed 
meaning, but <field name,field value> pairs on each line, separated by tabs. Each file can be 
loaded into Excel as a tab-delimited file and then filtered by record type, of first filtered 
by Unix 'grep'.

There is a field recType that can be used to filter the records by type, the values are:

* annotPred The record shows annotation-based majority voting prediction. If there were not 
enough hits to make a consensus prediction, assignment is left blank; if the hits did not 
provide the required majority vote even at a superkingdom level (in NCBI tree, superkindoms 
are Bacteria, Archaea and Eukaryots, we also call the Viruses a superkindom), the assignment 
is made to the root of the tree ("root").
* mgtPred The record shows MGTAXA prediction if those were provided as input, otherwise left 
blank.
* lcsPred The record shows the lowest common supernode (LCS) of MGTAXA-based and annotation-based 
predictions, if MGTAXA predictions were provided.
* annotRec The record shows a single BLASTP ORF assignment, one best hit per ORF, multiple records 
per contig.

For each record, a full main-rank lineage is reported, starting with bott_rank and bott_name 
fields that show the bottom-most rank and lineages assigned in a given record.

Optionally, a set of GFF3 and corresponding genomic diagrams in PDF and PNG formats are generated, 
in a subdirectory 'contigs", a set of files  per contig, named after the contig ID and the consensus 
annotation assignment.

Known issues
------------
------------

* Note on individual ORF annotations:

JCVI pipeline currently uses Uniref100 reference database for general-purpose BLASTP searches. 
100 means that sequences that are 100% pairwise identical over the length of a smaller sequence 
are clustered together. It is quite often that a cluster contains sequences from different 
organisms. In that case, Uniref100 reports the taxonomy of the lowest common supernode of 
those organisms, which is reflected in the ORF hit taxonomy produced by the annotation 
pipeline. This approach is reasonable because at 100% identity, there is no data to pick 
a leaf taxonomy of one sequence over another. An alternative approach that would pick the 
taxonomy according to a majority of leaf taxonomies of sequences in a cluster would be 
heavily biased by the over-representation of certain clades over others in the database. 
Examples of such clusters are http://www.uniprot.org/uniref/UniRef100_A7ZQU1 (assigned 
to Enterobacteriaceae) and even http://www.uniprot.org/uniref/UniRef100_P0CE57 (assigned 
to root because of a single bacteriophage sequence in the cluster).

In the annotation pipeline, the "Reviewed" records are split from "HighConfidence" records
if some of the Uniref cluster records come from SwissProt. They do not necesserily have better
BLAST scores than the remaining HighConfidence records. It the future, it might make sense to
pull the BLAST scores from btab file and select the best ORF record among HighConfidence and
Reviewed records based on the score.
"""
    
    def __init__(self,idSampToPred,taxaTree,taxaLevels,outDir):
        self.idSampSeen = set()
        self.idSampToPred = idSampToPred
        self.taxaTree = taxaTree
        self.taxaLevels = taxaLevels
        self.outDir = outDir
        makedir(self.outDir)
        self.outDirContigs = pjoin(outDir,"contigs")
        makedir(self.outDirContigs)
        self.outGff = None
        self.recGff = None
        self.gffFiles = {}
        self.annotCsvOut = openCsv(pjoin(outDir,"annot.csv"),
                mode="w",
                dialect="excel-tab",
                lineterminator="\n")
        self.predMetrics = dict(lcs=defdict(lambda:defdict(int)))

    def switchOutGff(self,idSamp,taxid):
        taxaTree = self.taxaTree
        outGff = self.outGff
        if outGff:
            outGff.close()
        taxaName = taxaTree.getNode(taxid).name if taxaTree and taxid > 0 else "unassigned"
        outFileRt = pjoin(self.outDirContigs,strToFileName("%s.%s" % \
                (idSamp,
                taxaName.replace(" ","_").replace("/","_"))))
        gffFile = outFileRt + ".gff3"
        if os.path.exists(gffFile):
            outGff = open(gffFile,"a")
        else:
            outGff = open(gffFile,"w")
            outGff.write(str(GFF3Header()))
        self.outGff = outGff
        self.recGff = GFF3Record(seqid=idSamp)
        self.gffFiles[gffFile] = {"idSamp":idSamp,"outFileRt":outFileRt}

    def closeGff(self):
        if self.outGff:
            self.outGff.close()
            self.outGff = None

    def closeCsv(self):
        if self.annotCsvOut.csvClose:
            self.annotCsvOut.csvFileStream.close()

    def cameraAnnot(self,inpAnnot,predMinLenSamp,annotClassifier):
        idSampToPred = self.idSampToPred
        idSampSeen = self.idSampSeen
        lastIdSamp = ""
        lastIdQ = ""
        picked = {}
        idSampAnnot = []
        for rec in inpAnnot:
            id_q = rec.parsed["id_q"]
            orf_cont = parseOrfOnContigId(id_q)
            idSamp = orf_cont["id_cont"]
            if idSamp in idSampToPred:
                pred = idSampToPred[idSamp]
                if pred[2] >= predMinLenSamp:
                    if not lastIdSamp:
                        idSampSeen.add(idSamp)
                        lastIdSamp = idSamp
                    if lastIdQ and lastIdQ != id_q:
                        for recType,parsed in picked.items():
                            idSampAnnot.append(parsed)
                        picked = {}
                    parsed = rec.parse()
                    parsed.update(orf_cont)
                    recType = parsed["type"]
                    if recType == "UNIREF_BLASTP":
                        #in annotation pipeline output, records are sorted
                        #by (id_q,unsorted by confidence,sorted by e-value within confidence)
                        if not recType in picked \
                                or picked[recType]["conf_x"] < parsed["conf_x"]:
                            picked[recType] = parsed
                        elif recType in picked \
                                and picked[recType]["conf_x"] == parsed["conf_x"]:
                                    print "DEBUG: Confidence tie: \n%s\n%s\n" % \
                                            (picked[recType],parsed)
                    if lastIdSamp != idSamp:
                        assert idSamp not in idSampSeen,\
                                "Discontinuous record order detected for sample ID: %s" % (idSamp,)
                        idSampSeen.add(idSamp)
                        if idSampToPred is not None:
                            mgtPred = idSampToPred[lastIdSamp]
                            mgtPredTaxid = mgtPred[1]
                            lenSamp = mgtPred[2]
                        else:
                            mgtPredTaxid = None
                            lenSamp = 0
                        #output accumulated annotations for last sample
                        annotPred = self.annotPredictAndWrite(idSamp=lastIdSamp,
                                annotRecs=idSampAnnot,
                                annotClassifier=annotClassifier,
                                mgtPredTaxid=mgtPredTaxid,
                                lenSamp=lenSamp)
                        annotPredNode = annotPred["annotPred"]["predNode"]
                        if annotPredNode:
                            gffPredTaxid = annotPredNode.id
                        else:
                            gffPredTaxid = None
                        self.gffWrite(idSamp=lastIdSamp,
                                annotRecs=idSampAnnot,
                                predTaxid=gffPredTaxid,
                                lenSamp=lenSamp)
                        lastIdSamp = idSamp
                        idSampAnnot = []
                    lastIdQ = id_q
        inpAnnot.close()

    def gffWrite(self,idSamp,annotRecs,predTaxid,lenSamp):
        self.switchOutGff(idSamp=idSamp,taxid=predTaxid)
        recGff = self.recGff
        iFeat = 1
        for parsed in annotRecs:
            recType = parsed["type"]
            recGff.type = "protein_match"
            recGff.start = parsed["start_orf"]
            recGff.end = parsed["end_orf"]
            recGff.strand = '+' if parsed["strand_orf"] > 0 \
                    else '-' if parsed["strand_orf"] < 0 \
                    else '.'
            #no attributes should be automatically carried forward from the previous values:
            recGff.attribs = GFF3Attributes()
            ats = recGff.attribs
            ats["ID"] = "%s_%s" % (parsed["id_s"],iFeat) #we should make ID unique for GFF
            ats["Name"] = "C: %(conf_x).1f TAXA: %(tax_name)s DESCR: %(descr)s" % parsed
            if len(ats["Name"]) == 0:
                ats["Name"] = ats["ID"]
            self.outGff.write(str(recGff))
            iFeat += 1

    def annotPredictAndWrite(self,idSamp,annotRecs,annotClassifier,mgtPredTaxid,lenSamp):
        taxaTree = self.taxaTree
        taxaLevels = self.taxaLevels
        annotPred = annotClassifier.predict(annotRecs=annotRecs)
        try:
            mgtPredNode = taxaTree.getNode(mgtPredTaxid)
        except KeyError:
            mgtPredNode = None
        self.writeCsvRecs(idSamp=idSamp,lenSamp=lenSamp,
                mgtPredNode=mgtPredNode,annotPred=annotPred,
                annotRecs=annotRecs)
        return dict(annotPred=annotPred)

    def writeCsvRecs(self,idSamp,lenSamp,mgtPredNode,annotPred,annotRecs):
        taxaTree = self.taxaTree
        taxaLevels = self.taxaLevels
        predMetrics = self.predMetrics
        annotPredNode = annotPred["predNode"]
        
        if annotPredNode and mgtPredNode:
            lcsPredNode = mgtPredNode.lcsNode(annotPredNode)
            lcsRec = {"distAnnotToLcs": annotPredNode.getDepth()-lcsPredNode.getDepth(),
                    "levAnnotToLcs": taxaLevels.getLevelPos(lcsPredNode.linn_level) - \
                            taxaLevels.getLevelPos(annotPredNode.linn_level),
                    "idLevLcs": taxaLevels.getLevelId(lcsPredNode.linn_level)}
        else:
            lcsPredNode = None
            lcsRec = {}
        
        #accumulate prediction metrics based on LCS node data
        if annotPredNode:
            skip = False
            applySkipRules = True
            if applySkipRules:
                #skip Actinobacteria (class) that we know MGT always misclassifes as
                #Polynucleobacter
                if annotPredNode.isUnder(taxaTree.getNode(1760)):
                    skip = True
                #skip if Annot assignment is no better than
                #a kingdom level to Proks. This covers very non-specific Annot
                #assignments to Proks which I suspect might be actually Euks or Viruses 
                #and detected as such by MGTAXA, but are 
                #not detected by BLAST due to very uneven DB coverage for marine pikoplankton
                #across superkingdoms
                #@todo check that Euk picoplankton genomes are in the BLAST DB used here
                elif ( not annotPredNode.isUnder(taxaTree.getNode(virTaxid)) ) and \
                        ( taxaLevels.getLevelId(annotPredNode.linn_level) >= \
                        taxaLevels.getLevelId("superkingdom") ):
                    skip = True
            if not skip:
                for level in ("order","class","phylum","superkingdom"):
                    if taxaLevels.getLevelId(annotPredNode.linn_level) <= \
                        taxaLevels.getLevelId(level):
                        predMetrics["lcs"][level]["nTest"] += 1
                        if lcsPredNode:
                            predMetrics["lcs"][level]["nMatch"] += \
                                    ( (lcsRec["idLevLcs"] <= taxaLevels.getLevelId("order")) or \
                                    lcsRec["distAnnotToLcs"] == 0 )

        annotCsvOut = self.annotCsvOut.csvFile
        for (botNode,outRecType,recRest) in [(annotPredNode,"annotPred",annotPred["ruleVars"]),
                (mgtPredNode,"mgtPred",{}),
                (lcsPredNode,"lcsPred",lcsRec)] + \
                        [ (annotPred["annotNodes"][id(annotRec)],"annotRec",annotRec) \
                        for annotRec in annotRecs ]:
            csvRec = OrderedDict([("idSamp",idSamp),("lenSamp",lenSamp)])
            csvRec["recType"] = outRecType
            csvRec["bott_rank"] = botNode.linn_level if botNode else None
            csvRec["bott_name"] = botNode.name if botNode else None
            levNames = taxaLevels.getLevelNames("ascend")
            if botNode:
                linFixed = taxaLevels.lineageFixedList(botNode,format="node")
            else:
                linFixed = [None]*len(levNames)
            for lev,node in zip(levNames,linFixed):
                csvRec[lev] = node.name if node else None
            for item in sorted(recRest.items()):
                csvRec[item[0]] = item[1]
            outRec = []
            for x in csvRec.items():
                outRec += [x[0],x[1]]
            annotCsvOut.writerow(outRec)

    def finPredMetrics(self):
        predMetrics = self.predMetrics
        predMetrics["lcs"] = dict(predMetrics["lcs"])
        for sub in predMetrics["lcs"].values():
            if "nTest" in sub:
                sub["accu"] = float(sub["nMatch"])/sub["nTest"]
        return predMetrics

    def graphics(self,gdFormat="pdf"):
        """Generate GFF files and graphics based on protein annotation for viral contigs.
        @param gdFormat file type for genomic digrams, one of "pdf","png","svg" 
        (svg output in genometools (v.1.3.5) incorrectly clips a lot on the right side of the diagram)"""
        for (gffFile,data) in self.gffFiles.items():
            grFile = "%s.%s" % (data["outFileRt"],gdFormat)
            pred = self.idSampToPred[data["idSamp"]]
            gr = GFF3Graphics(outFormat=gdFormat,width=max(pred[2]/10,800))
            try:
                gr(gffFile,grFile)
            except CalledProcessError, msg:
                print "Creating genome diagram from GFF3 file failed: %s" % (msg,)

    def writeDocs(self):
        strToFile(self._annotTaxaFileClassifierDoc,pjoin(self.outDir,"README.txt"))

if __name__ == '__main__':
    inpAnnot,inpPep,inpContigs,predTaxa,predMinLenSamp,outDir = sys.argv[1:]
    if inpPep == "None":
        inpPep = None
    if inpContigs == "None":
        inpContigs = None
    if predTaxa == "None":
        predTaxa = None
    assert not os.path.exists(outDir),"I will not overwrite data in existing directory: %s" \
            % (outDir,)
    if predTaxa:
        predTaxa = loadObj(predTaxa)
        idSampToPred = indexTuples(it.izip(predTaxa.idSamp,predTaxa.predTaxid,predTaxa.lenSamp))
    else:
        lenSamp = fastaLengths(inpContigs)
        #todo finish building idSampToPred
        idSampToPred = indexTuples(it.izip(lenSamp["id"],it.repeat(rejTaxid),lenSamp["len"]))
    
    if inpPep:
        inpPep = MappedPepFastaReader(inpPep)
        inpAnnot = CameraAnnotReaderPepIds(inpAnnot,pepMap=inpPep)
    else:
        inpAnnot = CameraAnnotReader(inpAnnot)
        
    predMinLenSamp = int(predMinLenSamp)
    #taxaTree = None
    taxaTree = loadTaxaTree()
    taxaLevels = TaxaLevels(taxaTree=taxaTree)
    #gffBtab(inpAnnot,predTaxa,taxaTree,predMinLenSamp,outDir)
    annotClassifier = AnnotTaxaClassifier(taxaTree=taxaTree,taxaLevels=taxaLevels)
    gffWriter = AnnotTaxaFileClassifier(idSampToPred=idSampToPred,taxaTree=taxaTree,taxaLevels=taxaLevels,outDir=outDir)
    gffWriter.cameraAnnot(inpAnnot=inpAnnot,predMinLenSamp=predMinLenSamp,annotClassifier=annotClassifier)
    gffWriter.closeCsv()
    gffWriter.closeGff()
    gffWriter.writeDocs()
    strToFile(str(gffWriter.finPredMetrics()),"annot.pred_metr.txt")
    gffWriter.graphics(gdFormat="pdf")
    gffWriter.graphics(gdFormat="png")

