### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Applicaion for  building and using IMM-based taxonomic classifier in the spirit of Phymm"""

from MGT.ImmApp import *
from MGT.Taxa import *
from MGT.App import *
from MGT.DirStore import *
from MGT.SeqDbFasta import *

class ImmClassifierApp(App):
    """App-derived class for building and using IMM-based taxonomic classifier in the spirit of Phymm.
    The main difference with Phymm is that we build IMMs for higher-level clades by pulling sequence
    data for subnodes.
    This class can be mostly viewed as imposing a TaxaTree structure onto ImmApp."""

    batchDepModes = ("predict","train")

    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        optChoicesMode = ("predict","train")
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=optChoicesMode,
            dest="mode",default="predict",help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            make_option(None, "--db-seq",
            action="store", type="string",dest="seqDb",help="Input SeqDbFasta path"),
            make_option(None, "--db-imm",
            action="store", type="string",default="imm",dest="immDb",help="Path to a collection of IMMs"),
            make_option(None, "--imm-seq-ids",
            action="store", type="string",default="imm-seq-ids",dest="immIdToSeqIds",help="File that maps IMM IDs to lists of seq IDs"),
            make_option(None, "--imm-ids",
            action="store", type="string",dest="immIds",help="File with list of IMM IDs to use in scoring and prediction. Default is all"+\
                    " IMMs from --imm-seq-ids"),
            make_option(None, "--inp-seq",
            action="store", type="string",dest="inpSeq",help="File with input FASTA sequence for prediction"),
            make_option(None, "--max-seq-id-cnt",
            action="store", type="int",default=100,dest="maxSeqIdCnt",help="Maximum number of training SeqDB IDs to propagate up from "+\
                    "every child of a given node"),
            make_option(None, "--out-dir",
            action="store", type="string",default="results",dest="outDir",help="Directory name for output score files"),
            make_option(None, "--out-score-comb",
            action="store", type="string",dest="outScoreComb",help="Output file for combined raw scores. Default "+\
                    "is 'combined'+ImmApp.scoreSfx inside --out-dir"),
        ]
        return Struct(usage = klass.__doc__+"\n"+\
                "%prog [options]",option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.setIfUndef("immIds",opt.immIdToSeqIds)
        opt.setIfUndef("outScoreComb",pjoin(opt.outDir,"combined"+ImmApp.scoreSfx))
    
    def initWork(self,**kw):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.seqDb = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
   

    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "train":
            return self.train(**kw)
        elif opt.mode == "predict":
            return self.predict(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def getSeqDb(self):
        opt = self.opt
        if self.seqDb is None:
            self.seqDb = SeqDbFasta.open(opt.seqDb) #"r"
            return self.seqDb

    def mapSeqToTree(self):
        """Assign list of SeqDB IDs to corresponding taxonomy tree nodes.
        In the current SeqDB version, each ID is a unique taxonomy id, so
        a one-element list will be assigned.
        @post Attribute 'leafSeqDbIds' is assigned to EVERY node and contains a list of IDs, possibly empty.
        The empty list will be a reference to a single shared object, so it should be treated as immutable"""
        taxaTree = self.getTaxaTree()
        seqDb = self.getSeqDb()
        taxaList = seqDb.getTaxaList()
        emptyList = []
        taxaTree.setAttribute("leafSeqDbIds",emptyList,doCopy=False)
        for taxid in taxaList:
            taxaTree.getNode(taxid).leafSeqDbIds = [ taxid ]


    def pickSeqOnTree(self,maxSeqIdCnt):
        """Assign to each node of the taxonomy tree the list of SeqDB IDs in the subtree under that node.
        It picks a medium count of IDs from each child node and itself, clipped by maxSeqIdCnt.
        """
        taxaTree = self.getTaxaTree()

        def actor(node):
            if node.isLeaf() and len(node.leafSeqDbIds)>0:
                node.pickedSeqDbIds = node.leafSeqDbIds
            else:
                chSeqIds = [c.pickedSeqDbIds for c in node.getChildren() if hasattr(c,"pickedSeqDbIds")]
                chSeqIds.append(node.leafSeqDbIds)
                chSeqIds = [ x for x in chSeqIds if len(x) > 0 ]
                if len(chSeqIds) > 0:
                    targLen = int(n.median([ len(x) for x in chSeqIds ]))
                    targLen = min(maxSeqIdCnt,targLen)
                    chSeqIds = [ sampleBoundWOR(x,targLen) for x in chSeqIds ]
                    node.pickedSeqDbIds = sum(chSeqIds,[])
                    assert len(node.pickedSeqDbIds) > 0

        taxaTree.visitDepthBottom(actor)

    def defineImms(self):
        taxaTree = self.getTaxaTree()
        immIdToSeqIds = {}
        for node in taxaTree.iterDepthTop():
            if hasattr(node,"pickedSeqDbIds"):
                immIdToSeqIds[node.id] = node.pickedSeqDbIds
        dumpObj(immIdToSeqIds,self.opt.immIdToSeqIds)

    def setupTraining(self):
        opt = self.opt
        self.mapSeqToTree()
        self.pickSeqOnTree(opt.maxSeqIdCnt)
        self.defineImms()

    def train(self,**kw):
        """Train all IMMs.
        Parameters are taken from self.opt
        """
        opt = self.opt

        self.setupTraining()

        optI = copy(opt)
        optI.mode = "train"

        imm = ImmApp(opt=optI)
        return imm.run(**kw)


    def predict(self,**kw):
        """Score with all IMMs and predict the taxonomy.
        Parameters are taken from self.opt
        @param inpSeq Name of the input multi-FASTA file to score
        @param outDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        optI = copy(opt)
        optI.mode = "score"
        imm = ImmApp(opt=optI)
        jobs = imm.run(**kw)
        #TMP:
        return jobs
        ##@todo Default method to make predictions based on scores (maybe)

        rOpt = copy(opt)
        rOpt.mode = "reduce-scores"
        rApp = ImmApp(opt=rOpt)
        kw = kw.copy()
        kw["depend"] = jobs

        return rApp.run(**kw)

    def reduceScores(self,**kw):
        """Reduce a matrix of combined scores to best positions on the taxonomic tree.
        Parameters are taken from self.opt
        @param outScoreComb name for file with combined scores
        @param outTaxaPred name for output file with predicted taxonomy
        """
        opt = self.opt
        sc = loadObj(opt.outScoreComb)

        #self.idImms = idImms
        #self.idScores = idScores
        #self.scores = scores
        raise NotImplementedError()
        pass
        
        idScores = None
        for (iImm,immId) in enumerate(opt.immIds):
            inpScoreFile = pjoin(opt.outDir,"%s%s" % (immId,self.scoreSfx))
            score = loadObj(inpScoreFile)
            if scores is None:
                scores = n.resize(score["score"],(len(score),len(opt.immIds)))
                idScores = score["id"].copy()
            else:
                assert n.all(idScores == score["id"])
            scores[:,iImm] = score["score"]
        dumpObj(ImmScores(idImms=opt.immIds,idScores=idScores,scores=scores),opt.outScoreComb)
    

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ImmClassifierApp)

