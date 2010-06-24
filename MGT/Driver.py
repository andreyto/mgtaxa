#!/usr/bin/env python
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


import sys, os
#sys.path.append(os.environ['CA_BINDIR'])

from MGT.CollectTaxa import *
from MGT.TreeSamplerApp import *
from MGT.ParamScanTestApp import *
from MGT.SeqImportApp import *
from MGT.Svm import *


def makeOpt():
    opt = Struct()
    opt.sampLen = 5000
    opt.cwd = pjoin(options.dataDir,"samp_%s" % opt.sampLen)
    opt.runMode = "inproc" #"inproc" #"batchDep"
    opt.rank = "genus"
    ## Subtrees starting at these nodes are considered "testing terminal" -
    ## entire subtree is always selected for testing. CV splits are also
    ## created along these boundaries.
    ## We use "sub..." ranks here because there are subspecies that do not have
    ## a species or genus node in their lineage.
    #opt.termTopRanks = ("species","subspecies") #("genus","species","subspecies","subgenus")
    opt.termTopRanks = (TreeSamplerApp.termTopRankLeaf,)
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

def makeOptParamScanTest():
    opt = makeOpt()
    opt.featName = "kmer_8" #"kmer_8a" #"kmerlad_4f_nnorm_rndGenGcFixed" #"kmerlad_4f" #"wd_2_10"
    clOpt = opt.setdefault("clOpt",Struct())
    psOpt = opt.setdefault("psOpt",Struct())
    tsOpt = opt.setdefault("tsOpt",Struct())
    ftOpt = opt.setdefault("ftOpt",Struct())
    #TMP:
    #clOpt.balanceTrainCounts = 100 #DBG
    clOpt.thresh = [-2000] #n.arange(-2.,1.01,0.1)
    clOpt.C = 1.
    clOpt.smlrConvergenceTol = 0.01 #faster, lower accuracy
    clOpt.smlrLm = 0.1
    clOpt.numTrainJobs = 100
    clOpt.MEM = 4000
    clOpt.saveConfMatr = True
    clOpt.exportConfMatr = True
    clOpt.method = "svm" #"svm" "knn" "smlr"
    clOpt.kernel = "lin" #"rbf"
    clOpt.knnK = 10
    ClassifierApp.fillWithDefaultOptions(clOpt)
    tsOpt.clOpt = deepcopy(clOpt)
    tsOpt.clOpt.MEM = 6000
    tsOpt.clOpt.numTrainJobs = 100
    tsOpt.runName = "%s-1" % clOpt.method
    pgen = ParamGridGen()
    #psOpt.params = pgen.add("C",[1]).grid() #(-10,10,2)(0,1,2)
    #psOpt.params = pgen.add("C",pgen.p2(6,16,2)).add("rbfWidth",pgen.p2(-4,10,2)).grid()
    psOpt.params = pgen.add("C",pgen.p2(-8,10,2)).grid() #(-10,10,2)(0,1,2)
    #psOpt.params = pgen.add("smlrLm",pgen.p2(-10,1,2)).grid() #(-10,10,2)(0,1,2)
    #psOpt.params = pgen.add("knnK",pgen.lin(1,20,2)).grid()
    #psOpt.params = pgen.add("knnMaxDist",pgen.lin(0.01,0.1,0.01)).grid()
    psOpt.clOpt = deepcopy(clOpt)
    psOpt.clOpt.MEM = 4000
    psOpt.clOpt.numTrainJobs = 10
    psOpt.runName = "%s-1" % clOpt.method
    #ftOpt.inSeq = "samp.rndGenGcFixed"
    ftOpt.featType = "kmer" #"wdh" "kmerlad" #"kmer"
    if ftOpt.featType == "wdh":
        ftOpt.revCompl = "merge" #"addcol" #"forward"
    elif ftOpt.featType == "kmerlad":
        ftOpt.kmerLen = 9
        ftOpt.revCompl = "addcol" #"forward"
        ftOpt.norm = NORM_POLICY.EXPECT | NORM_POLICY.EU_ROW
    elif ftOpt.featType == "kmer":
        ftOpt.kmerLen = 8
        ftOpt.revCompl = "merge" #"forward"
        ftOpt.norm = NORM_POLICY.FREQ | NORM_POLICY.EU_ROW
    ftOpt.balance = -2 #-2 #do not even shuffle (somehow it takes a while) - we do it anyway when making idlabs
    prOpt = opt.setdefault("prOpt",Struct())
    prOpt.clOpt = deepcopy(clOpt)
    prOpt.clOpt.MEM = 6000
    prOpt.clOpt.thresh = [-2000.]
    prOpt.runName = "1"
    return opt


def makeOptSampler():
    opt = makeOpt()
    opt.db = None
    return opt

def run_Sampler():
    opt = makeOptSampler()
    #modes = ["shred","split","label"]
    #modes = ["split","label"] # "shred" "split"
    modes = ["label"] # "shred" "split"
    opt.runMode = "inproc" #"batchDep"
    opt.MEM = 4000
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = TreeSamplerApp(opt=opt)
        jobs = app.run(depend=jobs)
    #return

def run_ParamScanTest():
    opt = makeOptParamScanTest()
    opt.runMode = "batchDep" #"inproc" #"batchDep"
    opt.MEM = 2000
    opt.cwd = pjoin(opt.cwd,opt.rank)
    #modes = ["feat","idlabs"]
    #modes = ["idlabs"]
    #modes = ["feat","idlabs","parscan"]
    modes = ["parscan"]
    #modes = ["feat","idlabs","test"]
    #modes = ["idlabs","test"]
    #modes = ["test"]
    #modes = ["train"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = ParamScanTestApp(opt=opt)
        jobs = app.run(depend=jobs)

def makeOptSeqImportApp():
    opt = makeOpt()
    opt.runMode = "inproc" #"inproc" #"batchDep"
    opt.cwd = pjoin(opt.cwd,opt.rank,"pred","asm11_Mar06")
    opt.inSeq = ["/usr/local/projects/GOSII/jmiller/asms/asm11_Mar06/gosII.scf.fasta.gz"]
    opt.minSampLen = 5000
    opt.inFormat = "ca"
    return opt

def run_SeqImportApp():
    opt = makeOptSeqImportApp()
    app = SeqImportApp(opt=opt)
    jobs = app.run()


def run_Pred():
    opt = makeOptParamScanTest()
    opt.runMode = "batchDep" #"inproc" #"inproc" #"batchDep"
    opt.cwd = pjoin(opt.cwd,opt.rank)
    optSI = makeOptSeqImportApp()
    opt.sampStorePred = optSI.cwd
    opt.prOpt.clOpt.thresh = [0.]
    opt.prOpt.name = "2"
    #modes = ["featPred","predict"]
    modes = ["predict"]
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = ParamScanTestApp(opt=opt)
        jobs = app.run(depend=jobs)


if __name__ == "__main__":
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
    #run_Sampler()
    run_ParamScanTest()
    #run_SeqImportApp()
    #run_Pred()
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
