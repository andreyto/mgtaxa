from shogun.Kernel import * #GaussianKernel, WeightedDegreeStringKernel, WeightedCommWordStringKernel
from shogun.Features import *
from shogun.Classifier import *
from shogun.PreProc import SortWordString, SortUlongString
import pdb
from MGT.PredProcessor import *
import random

def do_batch_linadd (train,test,C):
    print 'LibSVM batch'
    abet = DNA # IUPAC_NUCLEIC_ACID
    feats_train=StringCharFeatures(abet)
    feats_train.set_string_features(train['feature'])
    feats_test=StringCharFeatures(abet)
    feats_test.set_string_features(test['feature'])
    degree=20 #20

    kernel=WeightedDegreeStringKernel(feats_train, feats_train, degree)

    #C=0.017
    epsilon=1e-5
    tube_epsilon=1e-2
    num_threads=4
    labels=Labels(train['label'])

    svm=LibSVM(C, kernel, labels)
    svm.set_epsilon(epsilon)
    svm.set_tube_epsilon(tube_epsilon)
    svm.parallel.set_num_threads(num_threads)
    svm.train()

    kernel.init(feats_train, feats_test)

    #print 'LibSVM Objective: %f num_sv: %d' % \
    #    (svm.get_objective(), svm.get_num_support_vectors())
    #svm.set_batch_computation_enabled(False)
    #svm.set_linadd_enabled(False)
    #svm.classify().get_labels()

    svm.set_batch_computation_enabled(True)
    svm.set_linadd_enabled(True)
    pred = svm.classify()
    labpred = pred.get_labels()
    labtest = test['label']
    print numpy.bincount(numpy.sign(labpred)==labtest)

    #pdb.set_trace()

def seqToStringFeatures(train,test):

    feats_train=StringCharFeatures(DNA)
    feats_train.set_string_features(train['feature'].tolist())

    feats_test=StringCharFeatures(DNA)
    feats_test.set_string_features(test['feature'].tolist())

    return feats_train, feats_test

def seqToWordFeatures(train,test,order=7,gap=10):

    reverse=True

    charfeat_train, charfeat_test = seqToStringFeatures(train,test)

    feats_train=StringWordFeatures(charfeat_train.get_alphabet())
    feats_train.obtain_from_char(charfeat_train, order-1, order, gap, reverse)
    preproc=SortWordString()
    preproc.init(feats_train)
    feats_train.add_preproc(preproc)
    feats_train.apply_preproc()

    feats_test=StringWordFeatures(charfeat_test.get_alphabet())
    feats_test.obtain_from_char(charfeat_test, order-1, order, gap, reverse)
    feats_test.add_preproc(preproc)
    feats_test.apply_preproc()
    #pdb.set_trace()
    return feats_train, feats_test

def trainTest(svm,kernel,feats_train,feats_test):
    #pdb.set_trace()
    kernel.init(feats_train,feats_train)
    from MGT.Shogun import makePlattProb
    pb = makePlattProb(svm,nSplits=5)
    kernel.init(feats_train,feats_train)
    print "Training..."
    svm.set_batch_computation_enabled(True)
    svm.set_linadd_enabled(True)
    svm.train()
    print "Testing..."
    kernel.init(feats_train, feats_test)

    #print 'LibSVM Objective: %f num_sv: %d' % \
    #    (svm.get_objective(), svm.get_num_support_vectors())
    #svm.set_batch_computation_enabled(False)
    #svm.set_linadd_enabled(False)
    #svm.classify().get_labels()

    svm.set_batch_computation_enabled(True)
    svm.set_linadd_enabled(True)
    labTs = svm.classify()
    labTsPb = pb.apply(labTs) 
    pdb.set_trace()
    return labTsPb.get_labels()

def makeSvm(C,kernel,svmClass,labels=None):
    
    epsilon=1e-5
    tube_epsilon=1e-2
    num_threads=4
    if labels is not None:
        svm=svmClass(C, kernel, labels)
    else:
        svm=svmClass(C, kernel)

    svm.set_epsilon(epsilon)
    svm.set_tube_epsilon(tube_epsilon)
    svm.parallel.set_num_threads(num_threads)
    return svm

def runSvmBin(dataTrain,dataTest,kernelClass,C,order,gap):

    if kernelClass in (WeightedCommWordStringKernel,CommWordStringKernel):

        feats_train,feats_test = seqToWordFeatures(dataTrain,dataTest,order,gap)

        use_sign=False
        normalization=FULL_NORMALIZATION

        kernel=kernelClass(
            feats_train, feats_train, use_sign, normalization)

    elif kernelClass in (WeightedDegreeStringKernel,WeightedDegreePositionStringKernel):

        feats_train,feats_test = seqToStringFeatures(dataTrain,dataTest)

        degree = 7

        shift = 20

        if kernelClass is WeightedDegreeStringKernel:

            kernel=kernelClass(
                feats_train, feats_train, degree)

        else:

            kernel=kernelClass(
                feats_train, feats_train, degree)
            kernel.set_shifts(numpy.ones(len(dataTrain['feature'][0]), dtype='i4')*shift)

    labDataTrain = ((dataTrain['label'] - 1.5)*2).round()

    lab_train = Labels(labDataTrain)

    svm = makeSvm(C,kernel,LibSVM,lab_train)

    lab_pred = trainTest(svm,kernel,feats_train,feats_test)

    labDataPred = lab_pred/2. + 1.5

    for cutoff in numpy.arange(0.,2.,0.1):
    #for cutoff in numpy.arange(0.9,1.4,0.1):
        labpred = select([lab_pred <= -cutoff, lab_pred >= cutoff], [1, 2], default=0)
        labtest = dataTest['label']
        print "cutoff = ",cutoff
        perfMetrics(labtest,labpred)

def runSvmOne(dataTrain,dataTest,kernelClass,C,order,gap):
    labOne = int(dataTrain['label'][0])
    assert labOne != 0 and numpy.all(dataTrain['label'].astype('i4') == labOne)
    feats_train,feats_test = seqToWordFeatures(dataTrain,dataTest,order,gap)

    use_sign=False
    normalization=FULL_NORMALIZATION

    kernel=kernelClass(
        feats_train, feats_train, use_sign, normalization)

    svm = makeSvm(C,kernel,LibSVMOneClass)

    lab_pred = trainTest(svm,kernel,feats_train,feats_test)
    print lab_pred[:10]
    return
    pdb.set_trace()

    for cutoff in numpy.arange(0.,1.,1.):
        labpred = select([lab_pred <= cutoff, lab_pred > cutoff], [0, labOne], default=0)
        labtest = dataTest['label'].astype('i4')
        print "cutoff = ",cutoff
        perfMetrics(labtest,labpred)

def predSvmOne(dataTrain,dataTest,kernelClass,C,order,gap):
    labOne = int(dataTrain['label'][0])
    assert labOne != 0 and numpy.all(dataTrain['label'].astype('i4') == labOne)
    feats_train,feats_test = seqToWordFeatures(dataTrain,dataTest,order,gap)

    use_sign=False
    normalization=FULL_NORMALIZATION

    kernel=kernelClass(
        feats_train, feats_train, use_sign, normalization)

    svm = makeSvm(C,kernel,LibSVMOneClass)

    lab_pred = trainTest(svm,kernel,feats_train,feats_test)
    cutoff = 0.
    labpred = numpy.select([lab_pred <= cutoff, lab_pred > cutoff], [0, labOne], default=0)
    return labpred

def runSvmMult(dataTrain,dataTest,kernelClass,C,order,gap):
    feats_train,feats_test = seqToWordFeatures(dataTrain,dataTest,order,gap)

    use_sign=False
    normalization=FULL_NORMALIZATION

    kernel=kernelClass(
        feats_train, feats_train, use_sign, normalization)

    labDataTrain = (dataTrain['label'] - 1.).round()

    lab_train = Labels(labDataTrain)

    svm = makeSvm(C,kernel,LibSVMMultiClass,lab_train)

    lab_pred = trainTest(svm,kernel,feats_train,feats_test)

    labpred = (lab_pred + 1).round().astype('i4')
    labtest = dataTest['label'].astype('i4')
    perfMetrics(labtest,labpred)


def runSvmMultOne(dataTrain,dataTest,kernelClass,C,order,gap):
    feats_train,feats_test = seqToWordFeatures(dataTrain,dataTest,order,gap)

    use_sign=False
    normalization=FULL_NORMALIZATION

    kernel=kernelClass(
        feats_train, feats_train, use_sign, normalization)

    labDataTrain = (dataTrain['label'] - 1.).round()

    lab_train = Labels(labDataTrain)

    svm = makeSvm(C,kernel,LibSVMMultiClass,lab_train)

    lab_pred = trainTest(svm,kernel,feats_train,feats_test)

    labpred = (lab_pred + 1).round().astype('i4')
    labtest = dataTest['label'].astype('i4')

    labTrain = dataTrain['label'].astype('i4')

    labCnt = numpy.bincount(labTrain)

    for lab in xrange(len(labCnt)):
        if labCnt[lab] > 0:
            labpredOne = predSvmOne(dataTrain[labTrain == lab],dataTest,kernelClass,C,order,gap)
            #pdb.set_trace()
            labpred[labpred == lab] *= (labpred == labpredOne)[labpred == lab]

    perfMetrics(labtest,labpred)

def weighted_degree_string ():
    print 'WeightedDegreeString'

    feats_train=StringCharFeatures(DNA)
    feats_train.set_string_features(fm_train_dna)
    feats_test=StringCharFeatures(DNA)
    feats_test.set_string_features(fm_test_dna)
    degree=7 #20

    kernel=DegreeStringKernel(feats_train, feats_train, degree)

    #weights=arange(1,degree+1,dtype=double)[::-1]/ \
    #    sum(arange(1,degree+1,dtype=double))
    #kernel.set_wd_weights(weights)

    km_train=kernel.get_kernel_matrix()
    kernel.init(feats_train, feats_test)
    km_test=kernel.get_kernel_matrix(); pdb.set_trace()

def libsvm_oneclass ():
    print 'LibSVMOneClass'
    fm_train_real = nrnd.standard_normal((5,20))
    fm_test_real = nrnd.standard_normal((5,20))
    fm_test_real[:,:10] += 10
    feats_train=RealFeatures(fm_train_real)
    feats_test=RealFeatures(fm_test_real)
    width=2.1
    kernel=GaussianKernel(feats_train, feats_train, width)

    C=0.017
    epsilon=1e-5
    tube_epsilon=1e-2
    num_threads=1

    svm=LibSVMOneClass(C, kernel)
    svm.set_epsilon(epsilon)
    svm.set_tube_epsilon(tube_epsilon)
    svm.parallel.set_num_threads(num_threads)
    svm.train()

    kernel.init(feats_train, feats_test)
    #pred = svm.classify()
    #pdb.set_trace()
    print svm.classify().get_labels()

def transDegen(seq):
    abet = 'ACGT'
    nAbet = len(abet)
    s = numpy.fromstring(seq,dtype='S1')
    for i in xrange(len(s)):
        if s[i] not in abet:
            #print i,s[i],"->",
            s[i] = abet[nrnd.randint(nAbet)]
            #print s[i]
    return s.tostring()

def balance(data):
    cnt = numpy.bincount(data['label'].astype('i4'))
    targCnt = cnt[cnt > 0].min()
    dataSel = []
    for lab in xrange(len(cnt)):
        if cnt[lab] > 0:
            dataLab = data[data['label'].astype('i4') == lab]
            dataLab = dataLab[random.sample(xrange(len(dataLab)),targCnt)]
            dataSel.append(dataLab)
    dataSel = numpy.concatenate(dataSel)
    return dataSel

#libsvm_oneclass()
#sys.exit(0)

from MGT.Common import *
from itertools import izip

tr = {}
ts = {}

trans = ['C']*256
for c in 'ACTG':
    trans[ord(c)] = c
trans = ''.join(trans)

data = [None,None]

for seqF,dataInd in (("train.svm",0),("test.svm",1)):

    #inpSeq = FastaReader(pandaVirF)
    inpSeq = open(seqF,'r')
    labs = []
    seqs = []
    for line in inpSeq:
        lab, seq = line.split(None,1)
        lab = int(lab)
        if lab in (1,2,3,4):
            #labs.append((float(lab)-1.5)*2)
            labs.extend([lab])
            if seq[-1] == '\n':
                seq = seq[:-1]
            #seq = seq.translate(trans)
            seq = transDegen(seq)
            seqSplits = 4
            if seqSplits >= 2:
                sampLen = 1000/seqSplits
                iBeg = nrnd.randint(seqSplits-2)
                if dataInd == 0:
                    s = seq[sampLen*iBeg:sampLen*(iBeg+1)]
                    assert len(s) == sampLen
                    seqs.extend([s])
                else:
                    s = seq[-sampLen:]
                    assert len(s) == sampLen 
                    seqs.extend([s])
            else:
                seqs.append(seq)
    inpSeq.close()
    label = numpy.asarray(labs,dtype='i4')
    label[label==3] = 2
    label[label==4] = 2
    label = label.astype('f8')
    feature = numpy.asarray(seqs,dtype='O')
    data[dataInd] = numpy.rec.fromarrays((label,feature),names='label,feature')

tr,ts = data
#tr = tr[tr['label'] > 0]
tr = balance(tr)
print "Sample length: %s" % len(tr['feature'][0])
trInd = random.sample(xrange(len(tr)),800)
tr = tr[trInd]
print "len(tr) = %s counts(tr) = %s" % (len(tr), numpy.bincount(tr['label'].astype(int)))
#ts = balance(ts)
#tsInd = random.sample(xrange(len(ts)),700)
#tsMask = numpy.zeros(len(ts),bool)
#tsMask[tsInd] = True
#tsMask = numpy.logical_not(tsMask)
#ts = ts[tsMask]
print "len(ts) = %s" % len(ts)
#for C in numpy.power(2.,numpy.arange(-20,20,4)):
#for C in numpy.power(2.,numpy.arange(10,18,1)):
for C in (8,): #wdps
#for C in (256,): #wcw
#for C in (0.2,):
    print "C = ", C
    #do_batch_linadd(train=tr,test=ts,C=C)
    #weighted_comm_word_string(train=tr,test=ts,C=C)
    runSvmBin(dataTrain=tr,dataTest=ts,kernelClass=WeightedCommWordStringKernel,C=C,order=7,gap=10)
    #runSvmMultOne(dataTrain=tr,dataTest=ts,kernelClass=WeightedCommWordStringKernel,C=C,order=7,gap=10)
    #runSvmMult(dataTrain=tr,dataTest=ts,kernelClass=WeightedCommWordStringKernel,C=C,order=7,gap=10)
    #runSvmOne(dataTrain=tr[tr['label']==1],dataTest=ts,kernelClass=WeightedCommWordStringKernel,C=C,order=7,gap=10)

