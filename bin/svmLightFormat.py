"""Generate IdLabel file from SvmLight formatted sparse sequence file."""

from MGT.Svm import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-feat",
        action="store", type="string",dest="inFeat"),
        make_option("-o", "--out-lab",
        action="store", type="string",dest="outLab"),
        make_option("-s", "--split",
        action="store", type="int",dest="split",default=0,help="assign this split number to all data, default is %default"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()

data = loadSeqs(opt.inFeat,preProc=loadSeqPreprocDropFeat,genMissingId=True)
idFile = featIdFileNameDef(opt.inFeat)
if not os.path.exists(idFile):
    svmSaveId(data["id"],idFile)
split = n.zeros(len(data),dtype="i4")
split[:] = opt.split
idrecs = n.rec.fromarrays([data["id"],data["label"],split],names="id,label,split")
idLab = IdLabels(records=idrecs)
idLab.save(opt.outLab)

