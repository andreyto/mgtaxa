from MGT.Taxa import *
from MGT import Kmers
from itertools import izip

class SvmSparseFeatureWriterTxt(Kmers.SvmSparseFeatureWriterTxt):

    def write(self,label,feature):
        Kmers.SvmSparseFeatureWriterTxt.write(self,label,feature['values'],feature['indices'])

class SvmFastaFeatureWriterTxt:
    
    def __init__(self,out):
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        
    def close(self):
        self.out.close()

    def write(self,label,feature):
        self.out.write(">%d\n" % label)
        self.out.write(feature.tostring())
        self.out.write("\n")
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmStringFeatureWriterTxt:
    
    def __init__(self,out):
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        
    def close(self):
        self.out.close()

    def write(self,label,feature):
        self.out.write("%d " % label)
        if isinstance(feature,str):
            self.out.write(feature)
        else:
            self.out.write(feature.tostring())
        self.out.write("\n")
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmDenseFeatureWriterTxt:
    
    def __init__(self,out):
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = None
        
    def close(self):
        self.out.close()

    def write(self,label,values):
        if self.formatStr is None:
            self.formatStr = ' '.join( ("%d:%%g" % ind for ind in xrange(1,len(values)+1)) ) + '\n'
        self.out.write("%d " % (label,))
        self.out.write(self.formatStr % tuple(values))
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmDenseFeatureWriterCsv:
    
    def __init__(self,out,writeHeader=True,nFeat=None):
        assert not (writeHeader and nFeat is None)
        self.nFeat = nFeat
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = None
        if writeHeader:
            out.write("label,"+",".join([ "f%i" % iFeat for iFeat in xrange(nFeat) ])+"\n")
        
    def close(self):
        self.out.close()

    def write(self,label,values):
        if self.formatStr is None:
            self.formatStr = ",%g"*len(values)+'\n'
        self.out.write("%d" % (label,))
        self.out.write(self.formatStr % tuple(values))
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SVMLibLinear:
    
    def __init__(self,workDir):
        self.workDir = workDir
        self.modelFileRel = "model.svm"
        makedir(workDir)
        self.binDir = os.environ['SVM_LIBLINEAR_BIN']
        self.binTrain = os.path.join(self.binDir,'train')
        self.binPredict = os.path.join(self.binDir,'predict')
        
    def trainCmd(self,trainFile):
        #cmd = [self.binTrain,"-c","4","-e","0.1","-s","3",trainFile,self.modelFileRel]
        cmd = ["(zcat %s; echo '#'; zcat %s) |" % (trainFile,trainFile),
                self.binTrain,"-c","0.01","-e","1","-s","3","-",self.modelFileRel]
        return cmd
        #run(cmd, cwd=self.workDir)

        

class SVMLib:
    def __init__(self):
        import svm
        import cross_validation
        self.svm = svm
        self.cross_validation = cross_validation
        self.param = self.svm.svm_parameter(svm_type = self.svm.C_SVC, 
                                            kernel_type = self.svm.LINEAR, #self.svm.RBF,
                                            C = 2, gamma = 2, 
                                            cache_size=2000)
    
    def train(self,inpPath,testMode=False,testRecNum=10000):
        modelFile = inpPath + '.svm'
        kmersInp = KmerBinReader(inpPath)
        data = kmersInp.readAll()
        if testMode:
            data = numpy.concatenate((data[:testRecNum/2],data[-testRecNum/2:-1]))
        kmersInp.close()
        labels = data['taxid'].astype(numpy.float64)
        vals = data['vals'].astype(numpy.float64)
        print "labels: ", labels.shape, labels.dtype
        print "vectors: ", vals.shape, vals.dtype
        action = 'crossVal'
        if action == 'crossVal':
            print "Starting cross-validation with %s samples of rank %s" % (len(labels),len(vals[0]))
            self.cross_validation.do_cross_validation(vals,labels,self.param,5)
            print "Finished cross-validation"
        elif action == 'train':
            prob = self.svm.svm_problem(labels,vals)
            print "Starting training with %s samples of rank %s" % (len(labels),len(vals[0]))
            model = self.svm.svm_model(prob, self.param)
            print "Finished training"
            model.save(modelFile)
        else:
            raise ValueError("Unknown 'action' value: "+action)

class SVMMulticlass:
    
    def __init__(self,workDir='/export/tmp'):
        self.workDir = workDir
        self.sampleFileRel = "samples.svm"
        self.modelFileRel = "model.svm"
        makedir(workDir)
        
    def train(self,inpFile,testMode=False,testRecNum=10000):
        samples = open(os.path.join(self.workDir,self.sampleFileRel),'w', buffering=1024*1024)
        iRec = 0
        freqs = KmerFreqs(inpFile)    
        for rec in freqs.readValues():
            rec.write(samples)
            #samples.write(rec.asStr())
            samples.write("\n")
            iRec = iRec + 1
            if testMode and iRec >= testRecNum:
                break
        freqs.close()
        #important to flush, otherwise training program does not see the last records:
        samples.close()
        print "Starting training with %s samples of rank %s" % (iRec,len(rec.vals))        
        run(["svm_multiclass_learn","-c","1","-m","2000",self.sampleFileRel,self.modelFileRel], cwd=self.workDir)
        print "Finished training"



