from MGT.BioUtil import *

inpGb = pjoin(options.testDataDir,"gbank/problem_features.gbff")

def test_extractFeatureLocIndex():

    for rec in SeqIO.parse(inpGb,"genbank"):
        for feat in rec.features:
            if feat.type == "CDS":
                (pos,strand) = extractFeatureLocIndex(feat,format="range")
                print pos,strand
                (pos,strand) = extractFeatureLocIndex(feat,format="ind")
                print pos,strand

def test_FeatureLocIndex():
    outFeat = open("tmp.feat.fna","w")
    outNonFeat = open("tmp.nfeat.fna","w")
    for iRec,rec in enumerate(SeqIO.parse(inpGb,"genbank")):
        flind = FeatureLocIndex()
        for feat in rec.features:
            if feat.type == "CDS":
                flind.extract(feat)
        flind.toFasta(rec=rec,outFeat=outFeat,outNonFeat=outNonFeat)
    outFeat.close()
    outNonFeat.close()

#test_extractFeatureLocIndex()
test_FeatureLocIndex()

