from MGT.FeatIOCommon import *

class SvmFeatureWriterPkl:

    compressFormat = "gzip"
    
    def __init__(self,out,outId=None):
        if not hasattr(outId,'write'):
            if outId is None:
                assert isinstance(out,str)
                outId = open(featIdFileNameDef(out),'w', buffering=1024*1024)
        self.outId = outId
        self.out = PickleWriter(out,compressFormat=self.compressFormat)
        self.nOut = 0
        
    def close(self):
        self.out.close()
        self.outId.close()

    def write(self,label,feature,id=None):
        self.out.write((label,feature))
        if id is None:
            id = label
        self.outId.write("%s\n" % id)
        self.nOut += 1

    def numRec(self):
        return self.nOut

def SvmFeatureReaderPkl(inp):
    return PickleReader(inp,compressFormat=SvmFeatureWriterPkl.compressFormat)

