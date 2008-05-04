#!/usr/bin/env python
import sys, os
#sys.path.append(os.environ['CA_BINDIR'])

from MGT.CollectTaxa import *
#from MGT.Svm import *

dbSql = createDbSql()
#
db = TaxaCollector(dbSql=dbSql)
#db.tmp_loadSeqHdr()
db.rebuild()
#db.delDuplicateGiFromSeq()
#db.selectSeqIds()
#db.loadGiTaxNumpy()
#db.loadSeq()
#db.loadTaxNames()
#db.loadTaxLevels()
#db.loadTaxNodes()
#db.loadTaxCategories()
#db.loadRefseqAcc()
#db.selectTaxSet()
#db.loadTaxCategories()
#db.loadTaxNodes()
#db.loadTaxNodesMem()
#db.makeBlastAlias()
#db.mergeSelWithSeq(skipSeq=False)
#db.loadGiTaxSql()
#db.loadGiTaxNumpy()
#db.loadSeqToHdf()
#db.indexHdfSeq()
sys.exit(0)

nrnd.seed(1)
taxaSampler = TaxaSampler()
#taxaSampler.extractViralKmers()
#taxaSampler.extractViralSequence()
#taxaSampler.extractNonViralSequence()
#taxaSampler.makeSampleKmers()
#taxaSampler.makeViralVsNonViralSvm()
#taxaSampler.makeViralSvm()                     
#taxaSampler.writeViralFamilyLabels()
#updateKmerBinHeaderVersion(taxaSampler.kmerPathNonVir)
#kmers = KmerBinReader(taxaSampler.kmerPathNonVir)
#mergeKmerBin([taxaSampler.kmerPathVir,taxaSampler.kmerPathNonVir],taxaSampler.kmerPathVirNonVir)
#kmers = KmerBinReader(taxaSampler.kmerPath)
#ind = kmers.readIndTaxid()
#print numpy.sum(ind==0), len(ind)
#oldcnt = kmers.readCounts(update=True)
#newcnt = numpy.bincount(ind)
#assert numpy.all(oldcnt == newcnt)
#sys.exit(0)
#taxaSampler.subSampleKmersNonViral()
#taxaSampler.writeViralVsNonViralTrainingSet()
#taxaSampler.assignNodeCounts()
#sys.exit(0)
#trainSetPath = taxaSampler.trainSetPathVirFamily
#testSetPath = taxaSampler.testSetPathVirFamily
#kmerBalancedSample(inpPath=trainSetPath,outPath='tmp',medianRatio=0.2)
#kmerBinToSvmTxt(trainSetPath,trainSetPath+'.txt')
predRoot="lake_need_june"
predSeqPath = predRoot+".fas.gz"
#taxaSampler.sampleTestByPredictLength(predSeqPath)
#taxaSampler.sampleTestByPredictLengthSvm(predSeqPath)
#kmerBinToSvmTxt(taxaSampler.testSimPathVirFamily,taxaSampler.testSvmSimPathVirFamily+'.vnv',indLabelTaxid=numpy.ones(1000000,int))
#taxaSampler.makePredictSvm(predRoot)
taxaSampler.predictViral(predSeqPath,predRoot+'.6mers.vnv.svm.out',predRoot+'.7mers.v_f-tr.svm.out')
#trainSetPath = taxaSampler.trainSetPathVirNonVir
#kmerBinToSvmTxt(trainSetPath,trainSetPath+'.txt')

#svm = SVMMulticlass()
#svm = SVMLib()
#svm.train(trainSetPath,testMode=True,testRecNum=50000)
