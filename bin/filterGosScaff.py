"""Load and pickle Doug's scaffold-level taxonomy assignment file."""
from MGT.Taxa import *
from MGT.Common import *
from MGT.Shogun.Util import *
from MGT.Svm import *

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-file",
        action="store", type="string",dest="inFile"),
        make_option("-o", "--out-file",
        action="store", type="string",dest="outFile"),
        make_option("-d", "--doug-file",
        action="store", type="string",dest="dougFile"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


def filterScaff(sca,idPref=''):
    # a given level might have several records for a given scaffold
    # as long as we also require percent of peptides > 50, we always
    # choose only one
    sca = sca[
            logicalAnd(sca['level']==0,
                       sca['per_pep']>=95,
                       sca['length']>=5000,
                       sca['num_pep']*1000>=sca['length'])
            ]
    idsSel = (virTaxid,bacTaxid,archTaxid)
    idMap = {}
    counts = {}
    for idSel in idsSel:
        idRec = [ str(x) for x in sca[sca['taxid']==idSel]['id'] ]
        counts[idSel] = len(idRec)
        idSelSym = [ '%s%s_%s' % (idPref,idSel,idR) for idR in idRec ]
        idMap.update(dict(zip(idRec,idSelSym)))
    print "Filtering for: %s" % (sorted(counts.items()),)
    return idMap


opt,args = getProgOptions()

data = loadSeqs(opt.inFile)

scaff = loadObj(opt.dougFile)

idMap = filterScaff(scaff)

writer = SvmStringFeatureWriterTxt(opt.outFile)
counts = {}
for rec in data:
    if rec['id'] in idMap:
        id = idMap[rec['id']]
        writer.write(rec['label'],rec['feature'],id)
        idSel = id.split('_',1)[0]
        try:
            counts[idSel] += 1
        except KeyError:
            counts[idSel] = 1
print "Selected %s records: %s" % (writer.numRec(),sorted(counts.items()))
writer.close()

