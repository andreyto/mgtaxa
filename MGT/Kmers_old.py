from MGT.Common import *


class KmerFreqRecord:
    
    def __init__(self,label,vals):
        self.label = label
        self.vals = vals
        
    def asStr(self):
        """Return string representation accepted by both LibSVM and SVMLight binaries.
        Example: 3 1:0.43 3:0.12 9284:0.2"""
        return "%d " % (self.label,) + ' '.join( ( "%s:%s" % (i,v) for (i,v) in zip(range(1,len(self.vals)+1),self.vals) ) )
    
    def write(self,out):
        """Return string representation accepted by both LibSVM and SVMLight binaries.
        Example: 3 1:0.43 3:0.12 9284:0.2"""
        out.write("%d " % (self.label,))
        vals = self.vals
        for i in xrange(len(vals)):
            out.write("%s:%s " % (i+1,vals[i]))


    def __str__(self):
        return self.asStr()


class KmerBinHeader:
    def __init__(self,rootPath,nRec=0,nVal=0,recDtype=None):
        self.rootPath=rootPath
        self.nRec=nRec
        self.nVal=nVal
        self.recDtype=recDtype
        self.valPath=rootPath+'.val.gz'
        
    def countsPath(self):
        return self.rootPath+'.cnt'

#alows us to load older pickle dumps, need to run updateKmerBinHeaderVersion nevertheless
KmerBinData = KmerBinHeader

def updateKmerBinHeaderVersion(rootPath):
    oldHdr = loadObj(rootPath+'.hdr')
    hdr = KmerBinHeader(rootPath=oldHdr.rootPath,nRec=oldHdr.nRec,nVal=oldHdr.nVal,recDtype=oldHdr.recDtype)
    dumpObj(hdr,rootPath+'.hdr')


class KmerIndex:
    
    def __init__(self,name,rootPath):
        self.name = name
        self.rootPath = rootPath
                
    def load(self):
        return loadObj(self.fileName())

    def hasData(self):
        return os.path.isfile(self.fileName())
    
    def fileName(self):
        return self.rootPath + '.ind.' + self.name
    
    def save(self,data):
        dumpObj(data,self.fileName())

    def create(self):
        pass

    
class IndTaxid(KmerIndex):

    def create(self):
        kmerReader = KmerBinReader(rootPath=self.rootPath)
        #TODO: make kmerReader save max taxid in the header, use that 
        #here instead of dynamically resizing
        taxids = numpy.zeros(kmerReader.numRec(),dtype=numpy.int32)
        iRec = 0
        for batch in kmerReader.readBatches():
            taxids[iRec:iRec+len(batch)] = batch['taxid']
            iRec += len(batch)
            print "Read %s k-mer records out of %s" % (kmerReader.posRec(),kmerReader.numRec())
        kmerReader.close()
        self.save(taxids)
        return taxids


class KmerCounts:
    
    def __init__(self,rootPath=None,kmerHeader=None):
        assert (rootPath is None) ^ (kmerHeader is None)
        self.hdr = kmerHeader
        if self.hdr is None:
            self.hdr = KmerBinHeader(rootPath=rootPath)
        
    def load(self,update=False):
        fileName = self.hdr.countsPath()
        if self.hasCounts() and not update:
            return loadObj(fileName)
        else:
            reader = KmerBinReader(rootPath=self.hdr.rootPath)
            counts = self.create(reader)
            reader.close()
            return counts

    def hasCounts(self):
        return os.path.isfile(self.hdr.countsPath())
    
    def save(self,counts):
        """'counts' must have a property that it is zero for every id not present in the corresponding data file"""
        dumpObj(counts,self.hdr.countsPath())

    def create(self, kmerReader):
        #TODO: make kmerReader save max taxid in the header, use that 
        #here instead of dynamically resizing
        counts = numpy.zeros(100000,dtype=numpy.int32)
        maxTaxid = 0
        for batch in kmerReader.readBatches():
            for taxid in batch['taxid']:
                try:
                    counts[taxid] += 1
                except IndexError:
                    counts.resize(taxid*2)
                    counts[taxid] += 1
                maxTaxid = max(maxTaxid,taxid)
            print "Read %s k-mer records out of %s" % (kmerReader.posRec(),kmerReader.numRec())
        counts.resize(maxTaxid+1)
        self.save(counts)
        return counts
        

class KmerBinReader:
    
    def __init__(self,rootPath):
        self.hdr = loadObj(rootPath+'.hdr')
        self.inp = None
        #current input stream position in records
        self.iRec = 0
                         
    def readBatches(self,batchSize=10000):
        self._openInp()
        recDtype = self.hdr.recDtype
        nRec = self.hdr.nRec
        while True:
            if batchSize > nRec - self.iRec:
                batchSize = nRec - self.iRec
            if batchSize == 0:
                break
            recs = numpy.fromfile(self.inp, dtype=recDtype, count=batchSize)
            self.iRec += batchSize
            yield recs

    def readAll(self):
        return numpy.fromfile(self.inp, dtype=self.hdr.recDtype, count=self.hdr.nRec)
    
    def posRec(self):
        return self.iRec
    
    def numRec(self):
        return self.hdr.nRec
    
    def close(self):
        if self.inp is not None:
            self.inp.close()

    def readCounts(self,update=False):
        return KmerCounts(kmerHeader=self.hdr).load(update=update)

    def readIndTaxid(self):
        indTaxid = IndTaxid(name='taxid',rootPath=self.hdr.rootPath)
        if indTaxid.hasData():
            return indTaxid.load()
        else:
            return indTaxid.create()


    def _openInp(self):
        #numpy.fromfile does not recognize result of gzip.open() as a file object
        if self.inp is None:
            self.inp = Popen(("gzip -cd %s" % (self.hdr.valPath,)).split(),env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        
            

class KmerBinWriter:
    
    def __init__(self,rootPath):
        self.rootPath = rootPath
        self.hdr = KmerBinHeader(rootPath=self.rootPath,nRec=0,nVal=0,recDtype=None)
        self.out = None
    
    def writeBatch(self,recs):
        if len(recs) > 0:
            if self.out is None:
                self.hdr.nVal = len(recs[0]['vals'])
                self.hdr.recDtype = recs.dtype
                self.out = Popen("gzip -6 > %s" % (self.hdr.valPath,), shell=True, env=os.environ, bufsize=2**16, stdin=PIPE, close_fds=True).stdin
            recs.tofile(self.out)
            self.hdr.nRec += len(recs)
        
    def posRec(self):
        return self.hdr.nRec
    
    def numRec(self):
        return self.posRec()
        
    def close(self):
        if self.out is not None:
            self.out.close()
        dumpObj(self.hdr,self.rootPath+'.hdr')
        
    def saveCounts(self,counts=None):
        """'counts' must have a property that it is zero for every id not present in the corresponding data file"""        
        kmerCounts = KmerCounts(kmerHeader=self.hdr)
        if counts is None:
            kmerCounts.load(update=True)
        else:
            kmerCounts.save(counts)

    def saveIndTaxid(self,data=None):
        indTaxid = IndTaxid(name='taxid',rootPath=self.hdr.rootPath)
        if data is None:
            return indTaxid.create()
        else:
            indTaxid.save(data)
            

def mergeKmerBin(rootPaths,outRootPath):
    writer = KmerBinWriter(outRootPath)
    hdrOut = writer.hdr
    inpPaths = []
    for rootPath in rootPaths:
        reader = KmerBinReader(rootPath)
        hdr = reader.hdr
        inpPaths.append(hdr.valPath)
        hdrOut.nRec += hdr.nRec
        hdrOut.nVal = hdr.nVal
        hdrOut.recDtype = hdr.recDtype
        reader.close()
    #run('gunzip -c ' + ' '.join(inpPaths) + ' | gzip -6 -c > ' + hdrOut.valPath, shell = True)
    writer.close()


def kmerRecStr(rec):
    return "%d " % (rec['taxid']) + ' '.join( ( "%d:%.4f" % (i,v) for (v,i) in itertools.izip(rec['vals'],itertools.count(1)) if v != 0 ) )

class KmerToSvmWriter:
    
    def __init__(self,out,sparse=True):
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.formStr = None
        self.nOut = 0
        self.sparse = sparse
        
    def close(self):
        self.out.close()

    def write(self,rec):
        if self.sparse:
            self.writeSparse(rec)
        else:
            self.writeDense(rec)
        
    def writeDense(self,rec):
        self.out.write("%d " % (rec['taxid'],))
        vals = rec['vals']
        if self.formStr is None:
            self.formStr = ' '.join( ( "%d:%%.4f" % (i,) for i in xrange(1,len(vals)+1) ) ) + '\n'
        self.out.write(self.formStr % tuple(vals))
        self.nOut += 1

    def writeSparse(self,rec):
        self.out.write("%d " % (rec['taxid'],))
        vals = rec['vals']
        items = whereItems(vals,vals!=0)
        items['ind'] += 1
        strItems = ' '.join( ( "%d:%.4f" % tuple(item) for item in items ) )
        self.out.write(strItems)
        self.out.write("\n")
        self.nOut += 1

    def writeBatch(self,recs):
        for rec in recs:
            self.write(rec)

    def numRec(self):
        return self.nOut

def kmerBinToSvmTxt(inpPath,outFile,sparse=True,indLabelTaxid=None):
    kmersInp = KmerBinReader(inpPath)
    kmersOut = KmerToSvmWriter(outFile,sparse=sparse)
    for batch in kmersInp.readBatches(batchSize=1000):
        for rec in batch:
            #out.write(kmerRecStr(rec))
            #out.write('\n')
            if indLabelTaxid is not None:
                rec['taxid'] = indLabelTaxid[rec['taxid']]
            kmersOut.write(rec)
        print "Read %s k-mer records out of %s" % (kmersInp.posRec(),kmersInp.numRec())
    kmersInp.close()
    kmersOut.close()

def kmerTxtToSvmTxt(inpPath,outFile,sparse=True):
    inp = openGzip(inpPath,'r')
    out = open(outFile,'w', buffering=1024*1024)
    formStr = None
    for line in inp:
        rec = line.split(' ',5)
        taxid = 0
        print rec[:-1]
        vals = numpy.fromstring(rec[5],dtype=numpy.float32,sep=' ')
        #out.write(kmerRecStr(rec))
        #out.write('\n')
        out.write("%d " % (taxid,))
        if sparse:
            items = whereItems(vals,vals!=0)
            items['ind'] += 1
            strItems = ' '.join( ( "%d:%.4f" % tuple(item) for item in items ) )
            out.write(strItems)
            out.write("\n")
        else:
            if formStr is None:
                formStr = ' '.join( ( "%d:%%.4f" % (i,) for i in xrange(1,len(vals)+1) ) ) + '\n'
            out.write(formStr % tuple(vals))
        #for i in xrange(len(vals)):
        #    out.write("%d:%.4f " % (i+1,vals[i]))
        #out.write('\n')
            
    inp.close()
    out.close()


def kmerSampleByMask(inpPath,outPath,mask):
    inp = KmerBinReader(rootPath=inpPath)
    out = KmerBinWriter(rootPath=outPath)
    outTaxid = numpy.zeros(numpy.sum(mask == True),dtype=numpy.int32)
    iRec = 0
    iRecOut = 0
    for batch in inp.readBatches():
        maskSel = mask[iRec:iRec+len(batch)]
        batchSel = batch[maskSel == True]
        out.writeBatch(batchSel)
        outTaxid[iRecOut:iRecOut+len(batchSel)] = batchSel['taxid']
        iRec += len(batch)
        iRecOut += len(batchSel)
    inp.close()
    out.close()
    IndTaxid(name='taxid',rootPath=outPath).save(outTaxid)
    assert len(outTaxid) == iRecOut
    return outTaxid
    

def kmerBalancedSample(inpPath,outPath,medianRatio):
    inp = KmerBinReader(rootPath=inpPath)
    out = KmerBinWriter(rootPath=outPath)
    countsInp = inp.readCounts()
    medClassSize = numpy.median(countsInp[countsInp>0])
    targClassSize = int(medClassSize * medianRatio)
    print "balancedSample: selecting <= %s samples per class" % (targClassSize)
    countsOut = countsInp.clip(min=0,max=targClassSize)
    countsSel = numpy.zeros_like(countsOut)
    for batch in inp.readBatches():
        maskSel = numpy.zeros(len(batch),bool)
        for taxid in batch['taxid']:
            if countsSel[taxid] < countsOut[taxid]:
                #TODO:
                raise Error("Broken here")
                maskSel[taxid] = 1
                countsSel[taxid] += 1
        batchSel = batch[maskSel]
        out.writeBatch(batchSel)
    out.close()
    out.saveCounts(counts=countsSel)
    inp.close()
    assert numpy.all(countsSel == countsOut)


def makeKmersTxt(inpFile,outFile,kmerSize,kmerWindow=-1):
    genKmersExe = os.environ["PHYLA_GEN_KMERS_BIN"]
    run("zcat %s | %s -k %s -l %s | gzip -c > %s" % (inpFile,genKmersExe,kmerSize,kmerWindow,outFile),shell=True)

class KmerTxtReader:
    #tax_id, start on sequence, end on sequence, k-mer length, n of kmers, remaining columns are kmer frequencies 
    #(for example, k=6 has 2080 remaining columns)
    
    def __init__(self,inpFile,debug = True, nRecDbg = 10**3):
        self.debug = debug
        self.nRecDbg = nRecDbg
        self.inpFile = inpFile
        self.inp = gzip.open(inpFile,'r')
        self.countsFile = os.path.basename(inpFile)+".cnt"
    
    def readValues_new(self):
        debug = self.debug
        nRecDbg = self.nRecDbg
        iRec = 0
        valRank = -1
        while True:
            try:
                if iRec > 0:
                    vals = numpy.fromfile(self.inp, dtype=kmerDtype, count=valRank+5, sep=' ')
                else:
                    line = self.inp.readline()
                    vals = numpy.fromstring(line, dtype=float, count=-1, sep=' ')
                label = int(vals[0])
                vals = vals[5:]
                newRank = len(vals)
                kmerDtype = float
            except Exception, msg:
                print msg, line[:80], vals
                raise
            if valRank < 0:
                valRank = newRank
            else:
                assert valRank == newRank, "Old rank: %s, new rank %s" % (valRank,newRank)
                valRank = newRank
            iRec += 1
            if debug and iRec % nRecDbg == 0:
                print "Read %s k-mer vectors" % (iRec,)
            yield KmerFreqRecord(label=label,vals=vals)
        
    def readValues(self):
        debug = self.debug
        nRecDbg = self.nRecDbg
        iRec = 0
        valRank = -1
        for line in self.inp:
            try:
                vals = numpy.fromstring(line, dtype=float, count=-1, sep=' ')
                label = int(vals[0])
                vals = vals[5:]
                #rec = line.split()
                #label = int(rec[0])
                #vals = numpy.fromstring(','.join(rec[5:]), dtype=float, count=-1, sep=',')
                #vals = array(rec[5:],float)
                #vals = [ float(x) for x in rec[5:] ]
                newRank = len(vals)
            except Exception, msg:
                print msg, line[:80], vals
                raise
            if valRank < 0:
                valRank = newRank
            else:
                assert valRank == newRank, "Old rank: %s, new rank %s" % (valRank,newRank)
                valRank = newRank
            iRec += 1
            if debug and iRec % nRecDbg == 0:
                print "Read %s k-mer vectors" % (iRec,)
            yield KmerFreqRecord(label=label,vals=vals)
        
    def convertToBinary(self,debug=True,nRecDbg=10000):
        inp = openGzip(self.inpFile,'r')
        #inp = sys.stdin
        fileOutRoot = os.path.basename(self.inpFile)
        if fileOutRoot[-3:] == '.gz':
            fileOutRoot = fileOutRoot[:-3]
        fileOut = fileOutRoot + '.val.gz'
        #out = open(fileOut,'w')
        out = openGzip(fileOut,'w')
        iRec = 0
        line = inp.next()
        vals = numpy.fromstring(line, dtype=numpy.float32, count=-1, sep=' ')
        recDtype = numpy.dtype([('taxid','int32',1),
                                ('vals','float32',len(vals)-5)]) 
        recs = numpy.zeros(1,dtype=recDtype)
        rec = recs[0]
        while True:
            #rec = numpy.rec.fromarrays([vals[0:1].astype(numpy.int32),vals[5:]], names='taxid, vals')  
            rec['taxid'] = vals[0]
            rec['vals'][:] = vals[5:]
            rec.tofile(out)
            iRec += 1
            if debug and iRec % nRecDbg == 0:
                print "Read %s k-mer vectors" % (iRec,)
            #if iRec > 50000:
            #    break
            try:
                line = inp.next()
            except StopIteration:
                break
            vals = numpy.fromstring(line, dtype=numpy.float32, count=-1, sep=' ')
            
        out.close()
        inp.close()
        hdr = KmerBinHeader(rootPath=fileOutRoot,nRec=iRec,nVal=len(vals),recDtype=rec.dtype)
        dumpObj(hdr,fileOutRoot+'.hdr')
        
            
        
    def countValues(self,save=True,load=True):
        """If both 'save' and 'load' are True (the default), this will cache the results in a file."""
        if load and os.path.isfile(self.countsFile):
            return loadObj(self.countsFile)
        debug = self.debug
        nRecDbg = self.nRecDbg
        iniTaxaCount = 10**6
        count = numpy.zeros(iniTaxaCount,int)
        iRec = 0
        for line in self.inp:
            taxid = int(line.split(' ',1)[0])
            try:
                count[taxid] += 1
            except IndexError:
                count.resize((taxid*2,))
                count[taxid] += 1
            iRec += 1
            if debug and iRec % nRecDbg == 0:
                print "Read %s k-mer vectors" % (iRec,)
        if save:
            dumpObj(count,self.countsFile)
        return count
    
    def __del__(self):
        self.close()
    
    def close(self):
        self.inp.close()

def kmersTxtToBin(inpFile):
    reader = KmerTxtReader(inpFile)
    reader.convertToBinary()
    reader.close()
