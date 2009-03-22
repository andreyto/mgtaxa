#!/usr/bin/env python
import sys, os
#sys.path.append(os.environ['CA_BINDIR'])

from MGT.CollectTaxa import *
from MGT.TreeSamplerApp import *
from MGT.ParamScanTestApp import *
from MGT.SeqImportApp import *
from MGT.Svm import *

#dbSql = createDbSql()
#
#txcol = TaxaCollector(dbSql=dbSql)
#txcol.tmp_loadSeqHdr()
#txcol.rebuild()
#txcol.loadTaxLevelsSql()
#sys.exit(0)
#txcol.delDuplicateGiFromSeq()
#txcol.selectSeqIds()
#txcol.loadGiTaxNumpy()
#txcol.loadSeq()
#txcol.loadTaxNames()
#txcol.loadTaxLevels()
#txcol.loadTaxNodes()
#txcol.loadTaxCategories()
#txcol.loadRefseqAcc()
#txcol.selectTaxSet()
#txcol.loadTaxCategories()
#txcol.loadTaxNodes()
#txcol.loadTaxNodesMem()
#txcol.makeBlastAlias()
#txcol.mergeSelWithSeq(skipSeq=False)
#txcol.loadGiTaxSql()
#txcol.loadGiTaxNumpy()
#txcol.loadSeqToHdf()
#txcol.indexHdfSeq()

def makeOpt():
    opt = Struct()
    opt.sampLen = 5000
    opt.cwd = os.path.abspath("samp_%s" % opt.sampLen)
    opt.runMode = "inproc" #"inproc" #"batchDep"
    opt.rank = "genus"
    ## Subtrees starting at these nodes are considered "testing terminal" -
    ## entire subtree is always selected for testing
    ## We need "sub..." ranks because there are subspecies that do not have
    ## a species or genus node in their lineage.
    opt.termTopRanks = ("species","subspecies") #("genus","species","subspecies","subgenus")
    ## cross-validation splits will be done along these boundaries
    opt.splitRank = "species"
    opt.maxSamplesTrain = 60000
    ## label for reject group
    opt.labRj = 0
    ## label for background group
    opt.labBg = 1
    ## labels for other classes start from this
    opt.labCl = 2
    ## testing samples are under this split
    opt.splitIdTest = 1
    ## this split can be only used in training for final prediction
    opt.splitIdTrainPred = 2
    ## this and above is for splits used in cross-val
    opt.splitIdTrainCv = 3
    return opt

def makeOptSampler():
    opt = makeOpt()
    opt.db = None
    return opt

def run_Sampler():
    opt = makeOptSampler()
    #modes = ["shred","split","label"]
    modes = ["split","label"] # "shred" "split"
    opt.runMode = "batchDep" #"inproc" #"batchDep"
    opt.MEM = 4000
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = TreeSamplerApp(opt=opt)
        jobs = app.run(depend=jobs)
    #return

def run_ParamScanTest():
    opt = makeOpt()
    opt.runMode = "batchDep" #"inproc" #"batchDep"
    opt.MEM = 6000
    opt.cwd = pjoin(opt.cwd,opt.rank)
    #modes = ["feat","idlabs"]
    #modes = ["idlabs"]
    #modes = ["feat","idlabs","parscan"]
    #modes = ["parscan"]
    modes = ["test"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = ParamScanTestApp(opt=opt)
        jobs = app.run(depend=jobs)

def makeOptSeqImportApp():
    opt = makeOpt()
    opt.runMode = "inproc" #"inproc" #"batchDep"
    opt.cwd = pjoin(opt.cwd,opt.rank,"pred","asm11_Mar06")
    opt.inSeq = ["/usr/local/projects/GOSII/jmiller/asms/asm11_Mar06/9-terminator/gosII.scf.fasta"]
    opt.minSampLen = 5000
    opt.inFormat = "ca"
    return opt

def run_SeqImportApp():
    opt = makeOptSeqImportApp()
    app = SeqImportApp(opt=opt)
    jobs = app.run()


def run_Pred():
    opt = makeOpt()
    opt.runMode = "inproc" #"inproc" #"batchDep"
    opt.cwd = pjoin(opt.cwd,opt.rank)
    optSI = makeOptSeqImportApp()
    opt.sampStorePred = optSI.cwd
    #modes = ["featPred","predict"]
    modes = ["predict"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = ParamScanTestApp(opt=opt)
        jobs = app.run(depend=jobs)

#run_Sampler()
#run_ParamScanTest()
#run_SeqImportApp()
run_Pred()
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
