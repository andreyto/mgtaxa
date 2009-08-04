### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *
from MGT.Svm import *
from MGT.CRISPR import *
from MGT.FastaIO import *
from MGT.Common import *
import pdb

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-file",
        action="store", type="string",dest="inFile"),
        make_option("-o", "--out-spacers",
        action="store", type="string",dest="outSpacers"),
        make_option("-s", "--in-seq",
        action="store", type="string",dest="inSeq"),
        make_option("-p", "--out-seq",
        action="store", type="string",dest="outSeq"),
        make_option("-a", "--out-arrays",
        action="store", type="string",dest="outArrays"),
        make_option("-k", "--out-seq-left",
        action="store", type="int",dest="outSeqLeft",default=1000),
        make_option("-m", "--out-seq-right",
        action="store", type="int",dest="outSeqRight",default=1000),
        make_option("-l", "--min-spacer-len",
        action="store", type="int",dest="minSpacerLen",default=20),
        make_option("-u", "--max-spacer-len",
        action="store", type="int",dest="maxSpacerLen",default=60),
        make_option("-t", "--trim-spacer-len",
        action="store", type="int",dest="trimSpacerLen",default=0),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args


opt,args = getProgOptions()


inp = openCompressed(opt.inFile,'r')

FLD_IND_ID = 0
FLD_IND_POS = 1
FLD_IND_SPLEN = 4
FLD_IND_REP = -2
FLD_IND_SP = -1

assert opt.minSpacerLen > 0

dtype=[("idseq","O"),("idarr","i4"),("pos","i8"),("repeat","O"),("spacer","O")]
data = []
idarr = 0
idseqLast = ""
for line in inp:
    parts = [ part.strip() for part in line.split('\t') ]
    spLen = int(parts[FLD_IND_SPLEN])
    spLen = int(parts[FLD_IND_SPLEN])
    idseq = parts[FLD_IND_ID].split()[0]
    rep = parts[FLD_IND_REP]
    sp = parts[FLD_IND_SP]
    if idseq != idseqLast:
        idarr = 0
        idseqLast = idseq
    data.append((idseq,idarr,int(parts[FLD_IND_POS])-1,rep,sp))
    if spLen < 0:
        idarr += 1
data = n.asarray(data,dtype=dtype)
inp.close()

arr = {}
for iRec in xrange(len(data)):
    rec = data[iRec]
    key = (rec['idseq'],rec['idarr'].item())
    try:
        arr[key].append(rec)
    except KeyError:
        arr[key] = [ rec ]

for key,val in arr.items():
    iRep = 0
    newRecs = []
    for rec in val:
        # min spacer length check is done by piler
        #if iRep < len(val) - 1 and len(rec['spacer']) < opt.minSpacerLen + 2*opt.trimSpacerLen:
        #    continue
        newRecs.append(rec)
        if len(rec['spacer']) > opt.maxSpacerLen:
            assert iRep < len(val) - 1
            break
        iRep += 1
    arr[key] = newRecs
    if len(newRecs) < 2:
        del arr[key]

for key,val in arr.items():
    if crisprHasSimilarSpacers(val):
        print "Found putative array with similar spacers, deleting:\n%s\n%s" % (key,val)
        del arr[key]
    elif not crisprLengthRatioOk(val):
        print "Found putative array that violates length ratios, deleting:\n%s\n%s" % (key,val)
        del arr[key]


arrpos = {}
for key,val in arr.items():
    arrpos[key] = Struct(fullRange=(val[0]['pos'].item(),val[-1]['pos'].item()+len(val[-1]['repeat'].replace("-",""))),
            arr=val)

arrBySeq = {}
for key in sorted(arrpos.keys()):
    val = arrpos[key]
    keySeq = key[0]
    valSeq = val
    try:
        arrBySeq[keySeq].append(valSeq)
    except KeyError:
        arrBySeq[keySeq] = [ valSeq ]

data = Struct(arrByKey=arrpos,arrBySeq=arrBySeq)
dumpObj(data,opt.outArrays)

inpSeq = FastaReader(opt.inSeq)
outSeq = openCompressed(opt.outSeq,'w')

arrSeq = {}

for rec in inpSeq.records():
    idSeq = rec.getSimpleId()
    if idSeq in arrBySeq:
        seqArrays = arrBySeq[idSeq]
        seq = n.fromstring(rec.sequence().lower(),dtype='S1')
        for arr in seqArrays:
            arrRange = arr.fullRange
            arrReps = arr.arr
            startSeq = max(arrRange[0] - opt.outSeqLeft,0)
            startArr = arrRange[0] - startSeq
            endSeq = min(arrRange[1] + opt.outSeqRight,len(seq))
            if arrRange[1] > len(seq):
                print "erroneus array range %s for seq len %s, skipping array" % (arrRange,len(seq))
                continue
            for rep in arrReps:
                rseq = rep['repeat'].replace('-','')
                x,y = rep['pos'],rep['pos']+len(rseq)
                seq[x:y] = n.fromstring(seq[x:y].tostring().upper(),dtype='S1')
            # should output startSeq instead of arrRange[0]
            # convert all coordinates to absolute on seq and pickle
            outSeq.write(">%s_%s %s:%s %s %s\n" % (idSeq,arrRange[0],arrRange[0],len(seq)-arrRange[1],startArr,len(arrReps)))
            outSeq.write(seq[startSeq:endSeq].tostring())
            outSeq.write("\n")
            
outSeq.close()

outSp = openCompressed(opt.outSpacers,'w')
for key,val in arrpos.items():
    for rec in val.arr[:-1]:
        sp = rec['spacer']
        if len(sp) >= opt.minSpacerLen + 2*opt.trimSpacerLen:
            sp = sp[opt.trimSpacerLen:len(sp)-opt.trimSpacerLen]
            outSp.write(">%s_%s\n%s\n" % (rec['idseq'],rec['pos'],sp))
outSp.close()

