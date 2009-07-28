"""Generate IdLabel file from SvmLight formatted sparse sequence file or convert features between SvmLight and binary formats."""

from MGT.Svm import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-feat",
        action="store", type="string",dest="inFeat"),
        make_option(None, "--in-feat-format",
        action="store", type="choice",choices=featIOFormats,dest="inFeatFormat",default="txt"),
        make_option(None, "--out-feat",
        action="store", type="string",dest="outFeat",default=None),
        make_option(None, "--out-feat-format",
        action="store", type="choice",choices=featIOFormats,dest="outFeatFormat",default="pkl"),
        make_option("-o", "--out-lab",
        action="store", type="string",dest="outLab",default=None),
        make_option("-s", "--split",
        action="store", type="int",dest="split",default=0,help="assign this split number to all data, default is %default"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

if opt.outFeat is None:
    data = loadSeqs(opt.inFeat,preProc=loadSeqPreprocDropFeat,genMissingId=True,format=opt.inFeatFormat)
else:
    data = loadSparseSeqs(opt.inFeat,genMissingId=True,format=opt.inFeatFormat)
idFile = featIdFileNameDef(opt.inFeat)
if not os.path.exists(idFile):
    svmSaveId(data["id"],idFile)
if opt.outLab is not None:
    idrecs = IdLabels.makeRecords(len(data))
    split = n.zeros(len(data),dtype="i4")
    idrecs["split"] = opt.split
    idrecs["id"] = data["id"]
    idrecs["label"] = data["label"]
    idLab = IdLabels(records=idrecs)
    idLab.save(opt.outLab)
if opt.outFeat is not None:
    saveSparseSeqs(data,opt.outFeat,format=opt.outFeatFormat)

