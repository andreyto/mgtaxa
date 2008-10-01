"""Load and pickle Doug's scaffold-level taxonomy assignment file."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.Shogun.Util import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-o", "--out-file",
        action="store", type="string",dest="outFile"),
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

data = []
for o in inVir + [inAll,inNCBI]:
    d = loadSeqs(os.path.join(o.inDir,o.inFile))
    for rec in d:
        rec['id'] = o.inDir+'_'+rec['id']
    if o.inDir.startswith('NCBI'):
        d = balance(d,-1)
    data.append(d)

data = n.concatenate(data)

writer = SvmStringFeatureWriterTxt(opt.outFile)
for rec in data:
    writer.write(rec['label'],rec['feature'],rec['id'])
writer.close()

idPrefCnt = sorted(binCount([ tuple(id.split('_')[:-1]) for id in data['id'] ]).items())
print "Wrote %s records: %s" % (writer.numRec(),idPrefCnt)

