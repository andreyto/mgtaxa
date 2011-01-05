### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Applicaion for scaling ICM scores.
Scores produced by Glimmer ICMs have subtle dependencies on the length
of reference sequences used to train ICM models.
Testing with randomly generated references and queries (test/py/test_ImmAppInner.py)
demonstrated that scores have lower values for shorter references when
each reference was generated as a set of i.i.d. nucleotides. Inversly, when each
reference was generated as N copies of a randomly generated sequence, the
scores tended to have higher values for smaller N.
What this module does:
Given a collection of trained ICMs, it scores a collection of randomly generated queries
against each model, and uses the estimated mean and standard deviation to scale the score of
every real query during the classification phase.
The random queries are generated with a fixed length, and the statistics for the random
queries of the real query length are projected as mean and variance being proportional
to the length of the query sequence."""

from MGT.ImmApp import *
from MGT.Taxa import *
from MGT.App import *
from MGT.DirStore import *
from MGT.SeqDbFasta import *

def _randomSeq(length):
    # nseq samples, with 1 experiment in each returns nseq x 4 array 
    # with just one '1' in each row, 
    # then argmax along rows gives the winning nucleotide index
    # (seems there is a bug in numpy 1.3.0:
    # when length is int array element (nseq.dtype exists and is 'int64'),
    # nrnd.multinomial(1,p,length) returns a flat array (length,) instead of
    # (length,len(p))
    alph = n.fromstring("ATCG",dtype='S1')
    p = [0.25]*4
    return alph[nrnd.multinomial(1,p,int(length)).argmax(1)]

class ImmScalingApp(App):
    """App-derived class to generate random query sequences, score them and explore the distribution.
    """
    batchDepModes = ("score")

    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        optChoicesMode = ("generate","score")
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=optChoicesMode,
            dest="mode",help=("What to do, choice of %s, default is %%default" % (optChoicesMode,))),
            make_option(None, "--db-imm",
            action="store", type="string",default="imm",dest="immDb",help="Path to a collection of IMMs"),
            make_option(None, "--imm-ids",
            action="store", type="string",
            dest="immIds",help="File with list of IMM IDs to use in scoring and prediction. Default is all"+\
                    " in --db-imm"),
            make_option(None, "--pred-seq",
            action="store", type="string",
            dest="predSeq",help="FASTA file with random query sequences"),
            make_option(None, "--out-score-dir",
            action="store", type="string",dest="outScoreDir",help="Directory name for output score files"),
            make_option(None, "--out-score-comb",
            action="store", type="string",dest="outScoreComb",help="Output file for combined raw scores. Default "+\
                    "is 'combined'+ImmApp.scoreSfx inside --out-dir"),
            make_option(None, "--out-imm-scale-dir",
            action="store", type="string",default="scale",dest="outScaleDir",help="Directory name for output scaling files"),
            make_option(None, "--num-queries",
            action="store", type="int",default=10000,dest="numQueries",help="Number of random queries to generate and score"),
            make_option(None, "--query-len",
            action="store", type="int",default=5000,dest="queryLength",help="Length of a single query sequence"),
        ]
        return Struct(usage = klass.__doc__+"\n"+\
                "%prog [options]",option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.setIfUndef("cwd",os.getcwd())
        opt.setIfUndef("outScoreDir",pjoin(opt.cwd,"scores"))
        opt.setIfUndef("predSeq",pjoin(opt.cwd,"query.fna"))
        opt.setIfUndef("outScoreComb",pjoin(opt.outScoreDir,"combined"+ImmApp.scoreSfx))
    
    def initWork(self,**kw):
        opt = self.opt
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
   

    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "generate":
            return self.generate(**kw)
        elif opt.mode == "score":
            return self.score(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def generate(self,**kw):
        """Generate random query sequences.
        Parameters are taken from self.opt
        @param predSeq Output multi-FASTA file
        @param queryLength Length on generated sequences
        @param numQueries Number of generated sequences
        """
        opt = self.opt
        makeFilePath(opt.predSeq)
        writer = FastaWriter(opt.predSeq)
        lengths = n.ones(opt.numQueries,dtype=int)*opt.queryLength
        for iSeq,length in enumerate(lengths):
            seq = randomSeq(length)
            writer.record("%000i_%i" % (iSeq,length),seq)
        writer.close()

    def score(self,**kw):
        """Score with all IMMs and predict the taxonomy.
        Parameters are taken from self.opt
        @param predSeq Name of the input multi-FASTA file to score
        @param outScoreDir Directory name for output score files
        @param outScoreComb name for output file with combined scores
        """
        opt = self.opt
        optI = copy(opt)
        optI.inpSeq = opt.predSeq
        optI.outDir = opt.outScoreDir
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
    runAppAsScript(ImmScalingApp)

