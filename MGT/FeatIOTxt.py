### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.FeatIOCommon import *

from MGTX import kmersx
class SvmSparseFeatureWriterTxt(kmersx.SvmSparseFeatureWriterTxt):

    def __init__(self,out,outId=None):
        kmersx.SvmSparseFeatureWriterTxt.__init__(self,out)
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(featIdFileNameDef(out),'w', buffering=1024*1024)
        self.outId = outId

    def write(self,label,feature,id=None):
        kmersx.SvmSparseFeatureWriterTxt.write(self,int(label),feature['val'],feature['ind'])
        if id is None:
            id = label
        self.outId.write("%s\n" % id)

    def close(self):
        kmersx.SvmSparseFeatureWriterTxt.close(self)
        self.outId.close()

def loadSeqPreprocParseSparseOrig(lab,seq,id):
    feat = n.fromiter(( (int(ind),float(val)) 
        for ind,val in (entry.split(':') 
            for entry in seq.split()) ),dtype=MGTSparseRealFeatures.defDtype)
    return ([lab], [feat], [id])

_parseSparseTransTable = string.maketrans(':',' ')
#This one is about 70% faster than versions 1 & 2 below and 5x faster than the Orig
def loadSeqPreprocParseSparse(lab,seq,id):
    featDtype = MGTSparseRealFeatures.defDtype
    feat = n.fromstring(seq.translate(_parseSparseTransTable),sep=' ',dtype=featDtype["val"])
    feat = n.rec.fromarrays([feat[::2],feat[1::2]],dtype=featDtype)
    return ([lab], [feat], [id])

def loadSeqPreprocParseSparse1(lab,seq,id):
    feat = n.fromstring(seq.translate(_parseSparseTransTable),sep=' ',dtype='f8')\
            .view([("ind","f8"),("val","f8")]).astype([("ind","i4"),("val","f8")])
    return ([lab], [feat], [id])

def loadSeqPreprocParseSparse2(lab,seq,id):
    feat = n.rec.fromrecords([ entry.split(':') for entry in seq.split() ],dtype=[("ind","i4"),("val","f8")])
    return ([lab], [feat], [id])

class SvmFastaFeatureWriterTxt:
    
    def __init__(self,out,lineLen=None):
        if lineLen is None:
            lineLen = 1000
        if not hasattr(out,'write'):
            out = openCompressed(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.lineLen = lineLen
        
    def close(self):
        self.out.close()

    def write(self,label,feature,id=None):
        if label is None and id is not None:
            defLine = ">%s\n" % (id,)
        elif label is not None and id is None:
            defLine = ">%s\n" % (label,)
        else:
            defLine = ">%s|%s\n" % (id,label)
        self.out.write(defLine)
        if isinstance(feature,str):
            s = feature
        else:
            s = feature.tostring()
        lineLen = self.lineLen
        for x in range(0,len(s),lineLen):
            self.out.write(s[x:x+lineLen])
            self.out.write("\n")
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmStringFeatureWriterTxt:
    
    def __init__(self,out,outId=None):
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(featIdFileNameDef(out),'w', buffering=1024*1024)
        self.outId = outId
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        
    def close(self):
        self.out.close()
        self.outId.close()

    def write(self,label,feature,id=None):
        self.out.write("%s " % label)
        if isinstance(feature,str):
            self.out.write(feature)
        else:
            self.out.write(feature.tostring())
        self.out.write("\n")
        if id is None:
            id = label
        self.outId.write("%s\n" % id)
        self.nOut += 1

    def numRec(self):
        return self.nOut

def SvmStringFeatureReaderTxt(inp):
    if isinstance(inp,str):
        closeInp = True
        inp = openCompressed(inp,'rb')
    else:
        closeInp = False
    for line in inp:
        lab, seq = line.split(None,1)
        lab = float(lab)
        if seq[-1] == '\n':
            seq = seq[:-1]
        yield lab,n.fromstring(seq,dtype='S1').view(n.chararray)
    if closeInp:
        inp.close()

class SvmDenseFeatureWriterTxt:
    
    def __init__(self,out,outId=None):
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(featIdFileNameDef(out),'w', buffering=1024*1024)
        self.outId = outId
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = None
        
    def close(self):
        self.out.close()
        self.outId.close()

    def write(self,label,feature,id=None):
        if self.formatStr is None:
            self.formatStr = ' '.join( ("%d:%%g" % ind for ind in xrange(1,len(feature)+1)) ) + '\n'
        self.out.write("%s " % (label,))
        self.out.write(self.formatStr % tuple(feature))
        if id is None:
            id = label
        self.outId.write("%s\n" % id)
        self.nOut += 1

    def numRec(self):
        return self.nOut

class SvmDenseFeatureWriterCsv:
    
    def __init__(self,out,writeHeader=True,nFeat=None,sparseInput=False):
        assert not (writeHeader and nFeat is None)
        assert not (sparseInput and nFeat is None)
        self.sparseInput = sparseInput
        self.nFeat = nFeat
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = None
        if writeHeader:
            out.write("id,label,"+",".join([ "f%i" % iFeat for iFeat in xrange(nFeat) ])+"\n")
        
    def close(self):
        self.out.close()

    def write(self,label,feature,id=None):
        if self.formatStr is None:
            if self.sparseInput:
                nFeat = self.nFeat
            else:
                nFeat = len(feature)
                if self.nFeat is not None:
                    assert nFeat == self.nFeat
            self.formatStr = ",%g"*nFeat+'\n'
        if id is None:
            id = label
        self.out.write("%s,%d" % (id,label))
        if self.sparseInput:
            #sparse ind starts from 1
            x = n.zeros((self.nFeat,),dtype=feature.dtype["val"])
            x[(feature['ind']-1,)] = feature['val']
            feature = x
        self.out.write(self.formatStr % tuple(feature))
        self.nOut += 1

    def numRec(self):
        return self.nOut


class SvmDenseFeatureWriterXyz:
    """Output first three components of a dense real feature as atom coordinates of Tinker XYZ file.
    The output file can be loaded into some molecular viewer such as Pymol.
    The idea is to use fast 3D graphics and selection algebra of a molecular viewer
    in order to vizualize first three components of PCA or other dimensionality reduction procedure.
    Label is written as atom type."""
    
    def __init__(self,out,nSamp,outId=None):
        self.chemElem = "Ar"
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(featIdFileNameDef(out),'w', buffering=1024*1024)
        self.outId = outId
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        self.formatStr = "%6d  %-3s%12.6f%12.6f%12.6f%6d\n"
        #self.formatStr = "%i\t%s\t%f\t%f\t%f\%i\n"
        out.write("%6d  \n" % nSamp)
        
    def close(self):
        self.out.close()
        self.outId.close()

    def write(self,label,feature,id=None):
        num = self.nOut + 1
        if id is None:
            id = num
        self.out.write(self.formatStr % ((num,self.chemElem)+tuple(feature[:3])+(label,)))
        self.outId.write("%s\n" % id)
        self.nOut += 1

    def numRec(self):
        return self.nOut


