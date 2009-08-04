### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application interface to creation of features based on k-mer frequences out of biological sequences."""

from MGT.Shogun.Util import *
from shogun.Features import *
from MGT.Svm import *
from MGT.App import *

from MGT import Kmers
from MGT.Kmers import NORM_POLICY

__all__ = ["SeqFeaturesApp","NORM_POLICY"]


def otherVsRest(data,otherLabel):
    data['label'][:] = numpy.select([data['label'] == otherLabel],[1],default=2)
    return data



class SeqFeaturesApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", type="choice",choices=("default",),
            dest="mode",default="default"),
            make_option("-i", "--in-seq",
            action="store", type="string",dest="inSeq"),
            make_option(None, "--in-seq-format",
            action="store", type="choice",choices=featIOFormats,
            dest="inSeqFormat",default=defFeatIOFormat),
            make_option("-o", "--out-feat",
            action="store", type="string",dest="outFeat"),
            make_option(None, "--out-feat-format",
            action="store", type="choice",choices=featIOFormats,
            dest="outFeatFormat",default=defFeatIOFormat),
            make_option("-s", "--sigma",
            action="store", type="float",dest="sigma",default=200),
            make_option("-k", "--kmer-len",
            action="store", type="int",dest="kmerLen",default=2),
            make_option("-d", "--kmer-max-dist",
            action="store", type="int",dest="maxDist",default=10),
            make_option("-l", "--kmer-min-dist",
            action="store", type="int",dest="minDist",default=-1),
            make_option("-b", "--balance",
            action="store", type="int",dest="balance",default=-1),
            make_option("-u", "--other-group",
            action="store", type="int",dest="otherGroupLab",default=0),
            make_option("-f", "--feat-type",
            action="store", type="choice",choices=("wdh","kmer","kmerlad"),
            dest="featType",default="wdh"),
            make_option("-r", "--rev-compl",
            action="store", type="choice",choices=("merge","forward","addcol","addrow","reverse"),
            dest="revCompl",default="merge"),
            make_option(None, "--norm",
            action="store", type="string",dest="norm",default=None),
            make_option("-a", "--alphabet",
            action="store", type="choice",choices=("dna","protein"),
            dest="alphabet",default="dna"),
            make_option("-x", "--shred-len",
            action="store", type="int",dest="shredLen",default=-1),
            make_option("-y", "--shred-offset",
            action="store", type="int",dest="shredOffset",default=0),
            make_option("-z", "--shred-num",
            action="store", type="int",dest="shredNum",default=0),
        ]
        return Struct(usage = "Construct several types of word frequency features\n"+\
                "%prog [options]",option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        try:
            if options.norm is not None:
                options.norm = listOr(enumNamesToValues(options.norm,NORM_POLICY))
        except ValueError,AttributeError:
            parser.error("--norm value must be one or more of comma separated %s, received %s" % 
                    (enumNames(NORM_POLICY),options.norm))
    
    def doWork(self,**kw):
        opt = self.opt
        print "App options:\n", opt
        if opt.shredLen > 0:
            loadSeqPreproc = LoadSeqPreprocShred(sampLen=opt.shredLen,
                    sampNum=opt.shredNum,
                    sampOffset=opt.shredOffset)
        else:
            loadSeqPreproc = loadSeqPreprocIdent
        data = loadSeqs(opt.inSeq,preProc=loadSeqPreproc,format=opt.inSeqFormat)

        ##data = otherVsRest(data,2)
        print "Loaded " + showSvmDataCounts(data)
        if opt.balance >= -1:
            data = balance(data,opt.balance,labTargets={opt.otherGroupLab:-1})
        print "Balanced to " + showSvmDataCounts(data)
        ##TMP:
        ##data = splitStringFeat(data,750)
        if opt.alphabet == 'dna':
            shogAlpha = DNA
            applyToFeatData(data,transDegen)
        else:
            shogAlpha = PROTEIN
            assert opt.revCompl == 'forward'

        rcPolicyShog = WH_RC_FORWARD
        rcPolicyKmer = Kmers.RC_POLICY.DIRECT
        if opt.revCompl == "merge":
            rcPolicyShog=WH_RC_MERGE
            rcPolicyKmer=Kmers.RC_POLICY.MERGE
        elif opt.revCompl == "addcol":
            data = addRevComplCols(data)
        elif opt.revCompl == "forward":
            pass
        elif opt.revCompl == "addrow":
            data = addRevComplRows(data)
        elif opt.revCompl == "reverse":
            data = applyRevCompl(data)
        else:
            raise ValueError("Value %s for revCompl is not supported" % opt.revCompl)

        
        print "Computing features: " + showSvmDataCounts(data) 

        if opt.featType == "wdh":

            feat_char=StringCharFeatures(shogAlpha)
            data = applyToFeatData(data,lambda f: f.tostring())
            feat_char.set_string_features(data['feature'].tolist())

            feat_whd=WordHistogramFeatures()
            feat_whd.obtain_from_char(feat_char,opt.kmerLen,opt.maxDist,opt.sigma,opt.minDist,rcPolicyShog)
            print "WHD number of feature elements: %d" % feat_whd.get_num_elements()
            feat_sparse = feat_whd.get_sparse_real_features()
            print "WHD number feature dimensionality: %d" % feat_whd.get_num_features()

            lab = Labels(data['label'])
            feat_sparse.write_svmlight_file(opt.outFeat,lab)

            svmSaveId(data['id'],opt.outFeat+'.id')

        elif opt.featType in ("kmer","kmerlad"):
            maxSampLen = max( ( len(samp) for samp in data['feature'] ) )
            if opt.featType == "kmer":
                kmerCnt = Kmers.KmerSparseFeatures(sampLen=maxSampLen,
                        kmerLen=opt.kmerLen,
                        rcPolicy=rcPolicyKmer,
                        normPolicy=opt.norm)
            elif opt.featType == "kmerlad":
                kmerCnt = Kmers.KmerLadderSparseFeatures(sampLen=maxSampLen,
                        kmerLen=opt.kmerLen,
                        rcPolicy=rcPolicyKmer,
                        normPolicy=opt.norm)
            svmWriter = SvmSparseFeatureWriter(opt.outFeat,format=opt.outFeatFormat)
            for samp in data:
                feat = kmerCnt.kmerFrequencies(samp['feature'])
                svmWriter.write(int(samp['label']),feat,samp['id'])
            svmWriter.close()

        else:
            raise ValueError("Unknown feature type requested: %s" % opt.featType)

        print "Wrote " + showSvmDataCounts(data)


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(SeqFeaturesApp)

