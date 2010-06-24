### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Load and pickle Doug's scaffold-level taxonomy assignment file."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.Shogun.Util import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-o", "--out-file",
        action="store", type="string",dest="outFile"),
        make_option("-l", "--samp-len",
        action="store", type="int",dest="sampLen",default=5000),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


inVir = \
        [ Struct(inDir='GSIOVIR112',inFile='inp.samp'),
          Struct(inDir='GSIOVIR117',inFile='inp.samp'),
          Struct(inDir='GSIOVIR122',inFile='inp.samp') ]

inAll = Struct(inDir='GSIOMIC',inFile='inp.vm.samp')

inNCBI = Struct(inDir='NCBIVM',inFile='inp.samp.s.0')

opt,args = getProgOptions()

sampNumMIC = 1

_virSfx = "_%s" % virTaxid
def sampNumFuncMIC(label,seq,id):
    if id.endswith(_virSfx):
        return -1
    else:
        return sampNumMIC

# this will shred to all available samples for viral scaffolds, but only sampNumMIC max for others
sampPreProcMIC = LoadSeqPreprocShred(sampLen=opt.sampLen,sampNum=sampNumFuncMIC,sampOffset=0)
# this will shred to all available samples
sampPreProcOthers = LoadSeqPreprocShred(sampLen=opt.sampLen)

data = []
for o in inVir + [inAll,inNCBI]:
    if o.inDir.startswith('GSIOMIC'):
        sampPreProc = sampPreProcMIC
    else:
        sampPreProc = sampPreProcOthers
    d = loadSeqs(os.path.join(o.inDir,o.inFile),preProc=sampPreProc)
    for rec in d:
        rec['id'] = o.inDir+'_'+rec['id']
    if o.inDir.startswith('NCBI'):
        # balance by labels
        d = balance(d,1000)
    data.append(d)

data = n.concatenate(data)

writer = SvmStringFeatureWriterTxt(opt.outFile)
for rec in data:
    writer.write(rec['label'],rec['feature'],rec['id'])
writer.close()

idPrefCnt = sorted(binCount([ tuple(id.split('_')[:-1]) for id in data['id'] ]).items())
print "Wrote %s records: %s" % (writer.numRec(),idPrefCnt)

