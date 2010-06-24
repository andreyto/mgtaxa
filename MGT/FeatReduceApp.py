### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application interface to dimentionality reduction of real valued features."""

from MGT.Svm import *
from MGT.App import *
import mdp

__all__ = ["FeatReduceApp"]


class FeatReduceApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=("default",),
            dest="mode",default="default"),
            make_option("-p", "--method",
            action="store", type="choice",choices=("pca","kpca"),
            dest="method",default="pca"),
            make_option(None, "--kpca-estart",
            action="store", type="float",dest="kpcaEstart",default=2.),
            make_option(None, "--kpca-tau",
            action="store", type="float",dest="kpcaTau",default=None),
            make_option(None, "--kpca-repnum",
            action="store", type="int",dest="kpcaRepnum",default=10),
            make_option(None, "--kernel",
            action="store", type="choice",choices=("gauss","invdis","lin","pol","tanh"),
            dest="kernel",default="gauss"),
            make_option(None, "--invdis-eps",
            action="store", type="float",dest="invDisEps",default=1.),
            make_option(None, "--rbf-width",
            action="store", type="float",dest="rbfWidth",default=None),
            make_option(None, "--rbf-width-mult",
            action="store", type="float",dest="rbfWidthMult",default=1.),
            make_option("-i", "--in-feat",
            action="append", type="string",dest="inFeat"),
            make_option("-o", "--out-feat",
            action="store", type="string",dest="outFeat"),
            make_option("-a", "--labels",
            action="append", type="string",dest="labels"),
            make_option(None, "--pca-dim",
            action="store", type="float",dest="pcaDim",default=None),
            make_option(None, "--out-proc",
            action="store", type="string",dest="outProc",default=None),
            make_option(None, "--scale-before",
            action="store_true", dest="scaleBefore",default=False),
            make_option(None, "--scale-after",
            action="store_true", dest="scaleAfter",default=False),
            make_option(None, "--check-null-hyp",
            action="store", type="choice",choices=("varShuffle","elemRnd"),
            dest="checkNull",default=None),
            make_option("-b", "--balance",
            action="store", type="int",dest="balance",default=-2),
        ]
        return Struct(usage = "Reduce dimentionality of real valued features\n"+\
                "%prog [options]",option_list=option_list)

    def doWork(self,**kw):
        opt = self.opt
        print "App options:\n", opt
        if opt.method == "pca":
            return self.doPca(**kw)
        elif opt.method == "kpca":
            return self.doKPca(**kw)
        else:
            raise ValueError("Unknown opt.method value %s" % opt.method)

    def doPca(self,**kw):
        opt = self.opt
        x = self.prepFeatures()
        pcaNode = mdp.nodes.PCANode(output_dim = opt.pcaDim,svd=True)
        #pcaNode = mdp.nodes.NIPALSNode(output_dim = opt.pcaDim)
        #pcaNode = mdp.nodes.LLENode(k=50,output_dim = opt.pcaDim)
        #pcaNode = mdp.nodes.HLLENode(k=12,output_dim = opt.pcaDim)
        print "Starting Dimenstionality Reduction"
        pcaNode.train(x)
        pcaNode.stop_training()
        x = pcaNode.execute(x)
        print "Finished Dimenstionality Reduction"
        self.finalizeFeatures(x=x)
        #y = x
    
    def doKPca(self,**kw):
        from elefant.kernels import vector
        from elefant.kha.kha import KHA
        opt = self.opt
        x = self.prepFeatures()
        ## learning rate = estart * tau / (tau + iter number)
        ## there is a relation with width - if width gets larger,
        ## we need to increase the estart, or we get 1D strings out
        ## of data for any # of iterations - apparently, the iterative 
        ## procedure gets trapped in some local minimum
        ## estart = 10 worked for rbfWidthMult = 10, estart = 100 gave the nan
        estart = opt.kpcaEstart #2
        if opt.kpcaTau is None:
            tau = 3 * len(x)
        else:
            tau = opt.kpcaTau
        mu = 0.3
        eigen_type = 3

        repnum = opt.kpcaRepnum
        filename = None

        if opt.kernel == "gauss":
            kernel=vector.CGaussKernel()
            self.setRbfParam(kernel,x)
            kernel.RememberSquaredPart(x)
        elif opt.kernel == "invdis":
            kernel=vector.CInvDisKernel()
            kernel.epsilon = opt.invDisEps
            kernel.RememberSquaredPart(x)
        elif opt.kernel == "lin":
            kernel=vector.CLinearKernel()
        elif opt.kernel == "pol":
            kernel=vector.CPolynomialKernel()
            kernel.scale = 1.
            kernel.offset = 1.
            kernel.degree = 2
        elif opt.kernel == "tanh":
            kernel=vector.CTanhKernel()
            kernel.scale = 1.
            kernel.offset = 1.
        else:
            raise ValueError("Unknown value for opt.kernel: %s" % (opt.kernel,))
        print "Building kernel cache"
        kernel.RememberBasePart(x,x)

        rows = int(opt.pcaDim)
        assert n.allclose(rows,opt.pcaDim),"kpca support only integer opt.pcaDim parameters: %s" % opt.pcaDim

        print "Starting Dimenstionality Reduction"
        A, AK, errors = \
        KHA(T=x, repnum=repnum, rows=rows, eigen_type=eigen_type, estart=estart,
              tau=tau, mu=mu, lambd=0.99, verbose = False,
              kernel=kernel, filenamebase=filename) #kernel=kernel
        print "Finished Dimenstionality Reduction"
        self.finalizeFeatures(x=AK.T.A) #AK.T.A

    def loadLabels(self):
        self.idLab = loadIdLabelsMany(fileNames=self.opt.labels)
    
    def loadAsDense(self):
        sparse = loadSparseSeqsMany(self.opt.inFeat,idLab=self.idLab)
        return sparseToDenseSeqs(sparse)

    def saveAsSparse(self,x):
        saveDenseSeqsAsSparse(data=x,out=self.opt.outFeat)

    def setRbfParam(self,kernel,x):
        opt = self.opt
        if opt.rbfWidth is None:
            ##In Elefant, kernel = exp(-sigma * dot(x1,x2)), therefore scale is 1/(2*sigma^2)
            ##for a data set of points in
            ##an n-dimensional vector space drawn from a Gaussian distribution
            ##ep(-(||x1-x2||^2)/(2*sigma^2))
            #kernel.scale = float(1/x.std(0).mean())*0.1 #1e-4 #/ 128.0
            ## elefant method returns 1/median(||x1-x2||^2), so it must be robust to outliers
            kernel.EstimateSomeHyperparameters(x)
            rbfWidth = (0.5/kernel.scale)**0.5
        else:
            rbfWidth = opt.rbfWidth
        rbfWidth *= opt.rbfWidthMult
        kernel.scale = 0.5/rbfWidth**2
        print "kernel.scale = %s kernel.width = %s" % (kernel.scale, rbfWidth)
            

    def prepFeatures(self):
        opt = self.opt
        self.loadLabels()
        data = self.loadAsDense()
        print "Loaded " + showSvmDataCounts(data)
        if opt.balance >= -1:
            data = balance(data,opt.balance)
        print "Balanced to " + showSvmDataCounts(data)
        x = data["feature"]
        if opt.checkNull is not None:
            if opt.checkNull == "varShuffle":
                #shuffle each row seperately
                for rec in x:
                    nrnd.shuffle(rec)
            elif opt.checkNull == "elemRnd":
                x = nrnd.random_sample(x.shape)
        #TMP:
        #x = (x.T/n.sqrt((x*x).sum(1))).T
        x = x[:,:-4] #256+64+16+4; 340,84,20,4
        if options.debug > 0:
            print "Avg feature values[:10]: %s" % (x.mean(0)[:10],)
            print "Std feature values[:10]: %s" % (x.std(0)[:10],)
        if opt.scaleBefore:
            x = stdScaleAndCenter(x)
            if options.debug > 0:
                print "Avg feature values after scaling[:10]: %s" % (x.mean(0)[:10],)
                print "Std feature values after scaling[:10]: %s" % (x.std(0)[:10],)
        self.data = data
        return x

    def finalizeFeatures(self,x):
        opt = self.opt
        data = self.data
        if options.debug > 0:
            print "Avg feature values after PCA[:10]: %s" % (x.mean(0)[:10],)
            print "Std feature values after PCA[:10]: %s" % (x.std(0)[:10],)
        if opt.scaleAfter:
            x = stdScaleAndCenter(x)
            if options.debug > 0:
                print "Avg feature values after scaling[:10]: %s" % (x.mean(0)[:10],)
                print "Std feature values after scaling[:10]: %s" % (x.std(0)[:10],)
        data = makeDenseFeature(data["label"],x,data["id"])
        self.saveAsSparse(data)
        
    def normKmersByExpect(self,x):
        inp = openCompressed("/usr/local/projects/BIOINFO/WORK/data/db/ncbiMicrobial/kmer4.freq.gz")
        #@todo compute gc for each original sequence
        gc = []
        for line in inp:
            pass


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(FeatReduceApp)

