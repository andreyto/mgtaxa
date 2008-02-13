from MGT.Common import *

def dbClose(dbObj):
    print "Running 'atexit()' handler"
    dbObj.close()

class DbSQL:
    
    def __init__(self):
        atexit.register(dbClose, dbObj=self)

    def ddlIgnoreErr(self,*l,**kw):
        curs = self.cursor()
        try:
            curs.execute(*l,**kw)
        except StandardError, msg:
            print msg
            pass
        curs.close()
        
    def ddl(self,*l,**kw):
        curs = self.cursor()
        curs.execute(*l,**kw)
        curs.close()
        
    def execute(self,*l,**kw):
        curs = self.cursor()
        curs.execute(*l,**kw)
        return curs

    def executemany(self,*l,**kw):
        curs = self.cursor()
        curs.executemany(*l,**kw)
        curs.close()
 
    def dropTable(self,name):
        self.ddl("drop table if exists " + name)

    def dropIndex(self,name,table):
        self.ddl("drop index  if exists %s on %s" % (name,table))

    def connection(self):
        return self.con
    
    def cursor(self):
        return self.con.cursor()
    
    def commit(self):
        if hasattr(self,'con'):
            self.con.commit()

    def close(self):
        pass
   
class DbSQLLite(DbSQL):
    
    from pysqlite2 import dbapi2 as dbmod
    
    def __init__(self,dbpath,strType=str,dryRun=False):
        DbSQL.__init__(self)
        self.strType = strType
        self.dbpath = dbpath
        self.dryRun = dryRun
        self.con = self.dbmod.connect(self.dbpath)
        self.con.text_factory = self.strType


    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()


class DbSQLMy(DbSQL):
    
    import MySQLdb as dbmod
    
    def __init__(self,dryRun=False):
        DbSQL.__init__(self)        
        #self.dbmod.server_init(("phyla","--defaults-file=my.cnf"),('server',))
        self.con = self.dbmod.connect(unix_socket="/tmp/atovtchi.mysql.sock",
                                      host="localhost",
                                      db="phyla",
                                      user="root",
                                      passwd="OrangeN0",
                                      read_default_file="my.cnf",
                                      read_default_group="client")


    def close(self):
        self.commit()
        if hasattr(self,'con'):
            self.con.close()
        #self.dbmod.server_end()



def createDbSQL():
    #db = DbSQLLite("/export/atovtchi/test_seq.db")
    db = DbSQLMy()
    return db


class BulkInserter:
    
    def __init__(self,cursor,sql,bufLen):
        self.cursor = cursor
        self.bufLen = bufLen
        self.sql = sql
        self.buf = []
        
    def __call__(self,record):
        self.buf.append(record)
        if len(self.buf) == self.bufLen:
            self.flush()
            
    def flush(self):
        if len(self.buf) > 0:
            self.cursor.executemany(self.sql,self.buf)
            self.buf = []



class CountAggregateVisitor:
    
    def __call__(self,node):
        for child in node.children:
            node.data.kmerCnt += child.data.kmerCnt


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
    def __init__(self,rootPath,nRec,nVal,recDtype):
        self.rootPath=rootPath
        self.nRec=nRec
        self.nVal=nVal
        self.recDtype=recDtype
        self.valPath=rootPath+'.val.gz'

#compatibility with older pickle dumps
KmerBinData = KmerBinHeader

class KmerBinReader:
    
    def __init__(self,rootPath):
        self.hdr = loadObj(rootPath+'.hdr')
        #numpy.fromfile does not recognize result of gzip.open() as a file object
        self.inp = Popen(("gzip -cd %s" % (self.hdr.valPath,)).split(),env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        #current input stream position in records
        self.iRec = 0
                         
    
    def readValues(self,debug=True,nRecDbg = 10**3):
        print "DEBUG: ", self.valPath
        inp = gzip.open(self.valPath,'r')
        for iRec in xrange(self.nRec):
            rec = numpy.fromfile(inp, dtype=self.recDtype, count=1)
            if debug and iRec % nRecDbg == 0:
                print "Read %s k-mer vectors" % (iRec,)
            yield KmerFreqRecord(label=rec[0]['taxid'],vals=rec[0]['vals'])
        inp.close()
        
    def readBatches(self,batchSize=10000):
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
    
    def posRec(self):
        return self.iRec
    
    def numRec(self):
        return self.hdr.nRec
    
    def close(self):
        self.inp.close()
    
            

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
        self.out.close()
        dumpObj(self.hdr,self.rootPath+'.hdr')


class KmerTxtReader:
    
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
        inp = gzip.open(self.inpFile,'r')
        #inp = sys.stdin
        fileOutRoot = os.path.basename(self.inpFile)
        if fileOutRoot[-3:] == '.gz':
            fileOutRoot = fileOutRoot[:-3]
        fileOut = fileOutRoot + '.val.gz'
        #out = open(fileOut,'w')
        out = Popen("gzip -6 > %s" % (fileOut,),
                    shell=True, env=os.environ, bufsize=2**16, stdin=PIPE, close_fds=True).stdin
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


class TaxaSampler(PhyOptions):
    def __init__(self):
        PhyOptions.__init__(self)
        #kmerTxt = KmerTxtReader(self.kmerTestFile)
        #kmerTxt.convertToBinary()
        self.viralRootTaxid = 10239
        print "Loading taxonomy tree"
        #Pickling the taxa tree goes wrong apparently - memory consumption grows ~ 2x, and it takes longer than 
        #recreating from original flat file
        #self.taxaTree = objectDiskCacher(TaxaTree,os.path.basename(self.taxaNodesFile)+'.pkl')(ncbiDumpFile=self.taxaNodesFile)
        self.taxaTree = TaxaTree(ncbiDumpFile=self.taxaNodesFile,save=False,load=False)
        #kmerCounts = self.kmers.countValues(save=True,load=True)
        #print numpy.sum(kmerCounts)
        #return
        print "Finished loading taxonomy tree"
        mask = self.taxaTree.buildSubtreeMask(self.viralRootTaxid)
        kmers = KmerBinReader('6mers_3K')
        kmersSel = KmerBinWriter('6mers_3K.vir')
        for batch in kmers.readBatches(batchSize=10000):
            sel = batch[numpy.where(mask.take(batch['taxid']) > 0)]
            kmersSel.writeBatch(sel)
            print "Read %s k-mer records out of %s, selected %s" % (kmers.posRec(),kmers.numRec(),kmersSel.numRec())
        kmersSel.close()
        kmers.close()
        return
        #kmers = self.kmers
        for rec in kmers.readValues():
            pass
        #self.assignLeafCounts()        
        #self.taxaTree.writeLineage(sys.stdout)
        def node_printer(node):
            if node.data.rank == 'family':
                print node.data.kmerCnt, node.lineageRanksStr()
        #self.taxaTree.visitDepthTop(lin_printer,self.viralRootTaxid)
        viralRoot = self.taxaTree.getNode(self.viralRootTaxid)
        #rank_faker = RankFakerVisitor(topNode=viralRoot,lineage=viralRanksTemplate)
        #self.taxaTree.visitDepthTop(rank_faker,self.viralRootTaxid)
        #rank_reducer = RankReducerVisitor()
        #self.taxaTree.visitDepthTop(rank_reducer)
        count_aggregator = CountAggregateVisitor()
        self.taxaTree.visitDepthBottom(count_aggregator)
        self.taxaTree.visitDepthTop(node_printer,self.viralRootTaxid)
    
    def assignLeafCounts(self):
        kmerCounts = self.kmers.countValues(save=True,load=True)
        nodes = self.taxaTree.getNodesDict()
        for node in nodes.itervalues():
            try:
                x = kmerCounts[node.data.taxid]
            except IndexError:
                x = 0
            node.data.kmerCnt = x
        
        

class SVMLib:
    import svm
    def __init__(self):
        self.param = self.svm.svm_parameter(kernel_type = self.svm.RBF,C = 10,cache_size=2000)
    
    def train(self,inpFile,modelFile,testMode=False,testRecNum=10000):
    #tax_id, start on sequence, end on sequence, k-mer length, n of kmers, remaining columns are kmer frequencies 
    #(for example, k=6 has 2080 remaining columns)
        labels = []
        samples = []
        freqs = KmerFreqs(inpFile)
        for rec in freqs.readValues():
            labels.append(rec.label)
            samples.append(rec.vals)
            if testMode and len(labels) >= testRecNum:
                break
        freqs.close()
        prob = self.svm.svm_problem(labels,samples)
        print "Starting training with %s samples of rank %s" % (len(labels),len(samples[0]))
        model = self.svm.svm_model(prob, self.param)
        print "Finished training"
        model.save(modelFile)

class SVMMulticlass:
    
    def __init__(self,workDir='/export/tmp'):
        self.workDir = workDir
        self.sampleFileRel = "samples.svm"
        self.modelFileRel = "model.svm"
        makedir(workDir)
        
    def train(self,inpFile,modelFile,testMode=False,testRecNum=10000):
        samples = open(os.path.join(self.workDir,self.sampleFileRel),'w', buffering=1024*1024)
        iRec = 0
        freqs = KmerFreqs(inpFile)    
        for rec in freqs.readValues():
            rec.write(samples)
            #samples.write(rec.asStr())
            samples.write("\n")
            iRec = iRec + 1
            if testMode and iRec >= testRecNum:
                break
        freqs.close()
        #important to flush, otherwise training program does not see the last records:
        samples.close()
        print "Starting training with %s samples of rank %s" % (iRec,len(rec.vals))        
        run(["svm_multiclass_learn","-c","1","-m","2000",self.sampleFileRel,self.modelFileRel], cwd=self.workDir)
        print "Finished training"

class DbSeqSource(PhyOptions):
    """Database of sequence source for training the classiffier"""
    
    def __init__(self,dbSql):
        
        PhyOptions.__init__(self)
        
        self.dbSql = dbSql
        self.ncbiDbs =  (
                         Struct(id='g',db='refseq_genomic'),
                         Struct(id='o',db='other_genomic'),
                         Struct(id='n',db='nt'),
                         Struct(id='h',db='htgs'),
                         Struct(id='w',db='wgs')
                         )


    def makeBlastAlias(self,idList=None,dbNames=None):
        from glob import glob
        from datetime import datetime
        from textwrap import dedent
        if dbNames is None:
            dbNames = [ db.db for db in self.ncbiDbs ]
        directDbRefs = False
        #comment out GILIST field in alias file if idList is None
        if idList is None:
            comGILIST = '#'
        else:
            comGILIST = ''
        dbNameAlias = self.srcDbNameAlias
        idListBin = dbNameAlias+'.gil'
        aliasStr = dedent("""\
        #
        # Alias file created %s
        #
        #
        TITLE Custom Unified DB for Phyla
        #
        DBLIST %%s
        #
        %sGILIST %s
        #
        #OIDLIST
        #
        """ % (datetime.today().ctime(),comGILIST,idListBin))
        cwd = os.getcwd()
        try:
            os.chdir(self.blastDataDir)
            if directDbRefs:
                chunks = []
                for rec in dbNames:
                    chunks += sorted(list(set(
                                    ( y.group(1) for y in  
                                      (re.match('('+rec.db+'\.[0-9]+)\..*',x) 
                                       for x in glob(rec.db+'.*')) 
                                    if y is not None )
                                    )))
            else:
                chunks = dbNames
            strToFile(aliasStr % ' '.join(chunks),dbNameAlias+'.nal')
        finally:
            os.chdir(cwd)
        
    def mergeSelWithSeq(self,skipSeq=False):
        #from itertool import izip
        outFasta = gzip.open(self.selFastaFile,'w',compresslevel=4)
        inpDump = open(self.selDumpFile,'r')
        selGiFile = os.path.abspath(self.selGiFile)
        fldsDump = "gi,taxid,src_db,kind,project,cat,stage,src_type,iid,seq_len,divid,rank".split(',')
        pipe = Popen(("fastacmd -i %s -d %s" % (selGiFile,self.srcDbNameAlias)).split(), 
                     cwd=self.blastDataDir, env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        #inpFasta = readFastaRecords(pipe,readSeq=True)
        FGI_SKIP  = 0x01
        FGI_WRITE = 0x02
        FGI_MISM  = 0x04
        giSeen = {}
        iRec = 0
        skip = True
        mismatchRun = 0 # how many gi mismatches in a row
        for line in pipe:
            try:
                if line.startswith(">"):
                    skip = False
                    #header line can be:
                    #>gi|23455713|ref|NC_004301.1| Enterobacteria phage FI, complete genome >gi|15183|emb|X07489.1| Bacteriophage SP genomic RNA
                    headers = line.split('>')
                    gi2 = -1
                    if len(headers) > 2:
                        hdr2 = headers[2]
                        if hdr2.startswith('gi|'):
                            gi2 = int(hdr2.split('|',2)[1])
                            if not giSeen.has_key(gi2):
                                giSeen[gi2] = 0
                    (gifld,gi,accfld,acc,txt) = line[1:].split('|',4)
                    assert gifld == 'gi'
                    gi = int(gi)
                    valsDump = inpDump.readline().rstrip('\n').split(' ') #empty string fields will be ok
                    giDump = int(valsDump[0])
                    taxidDump = int(valsDump[1])
                    try:
                        lineage = self.taxaTree.getNode(taxidDump).lineageRanksStr()
                    except KeyError:
                        lineage = 'NULL'
                        print "Warning: Lineage not found for taxid %s" % (taxidDump,)
                    if giDump == gi:
                        line = ">gi|%s|%s|%s|" % (gi,accfld,acc) + \
                            ''.join(["%s:%s " % (fld,val) for (fld,val) in zip(fldsDump[1:],valsDump[1:])]) + \
                            "lineage:%s " % (lineage,) + \
                            txt
                        if gi2 > 0:
                            giSeen[gi2] |= FGI_WRITE
                        mismatchRun = 0                            
                    elif giDump == gi2:
                        giSeen[giDump] |= FGI_SKIP
                        skip = True
                        mismatchRun = 0
                    else:
                        if mismatchRun >= 10:
                            raise ValueError("Mismatch between FASTA and SQl DUMP input streams:" + \
                                "fastaTitle = %s valsDump = %s" % (line,' '.join(valsDump)))
                        else:
                            if not giSeen.has_key(giDump):
                                giSeen[giDump] = 0
                            giSeen[giDump] |= (FGI_SKIP | FGI_MISM)
                            skip = True
                            mismatchRun += 1
                            print "GI mismatch for ", giDump
                    if iRec % 10000 == 0:
                        print "Done %s records" % (iRec,)
                    iRec += 1
                elif skipSeq:
                    skip = True
                if not skip:
                    outFasta.write(line)
            except:
                print "Exception with input line: ", line
                pipe.close()
                raise
        pipe.close()
        inpDump.close()
        outFasta.close()
        print 'giSeen = \n', sorted(giSeen.items())
        print 'len(giSeen) = ',len(giSeen)
        for (gi,val) in sorted(giSeen.items()):
            if val & FGI_SKIP and not val & FGI_WRITE:
                print "%s never written %s" % (gi,val)
            if val & FGI_MISM:
                print "%s mistamtch %s" % (gi,val)


    def loadSeq(self):
        self.loadGiTaxPickled()
        #self.createTableSeq()
        #self.loadSeqNCBI(self.ncbiDbs[0])
        #self.loadSeqNCBI(self.ncbiDbs[1])
        #self.loadSeqNCBI(self.ncbiDbs[2])
        #self.loadSeqNCBI(self.ncbiDbs[3])
        #self.loadSeqNCBI(self.ncbiDbs[4])
        #self.dbSql.ddl("create index taxid on seq(taxid)")
    
    def createTableSeq(self):
        self.dbSql.dropIndex("taxid","seq")
        self.dbSql.dropTable("seq")
        self.dbSql.execute(
        """
        create table seq
        (
        iid integer auto_increment primary key ,
        gi bigint,
        taxid integer,
        src_db varchar(1),
        project varchar(4),
        seq_len bigint,
        acc varchar(20),
        kind varchar(2),
        seq_hdr varchar(%s)
        )
        engine myisam
        """ % (self.fastaHdrSqlLen,))
        
    def loadSeqNCBI(self,db):
        pipe = Popen(("fastacmd -D 1 -d " + db.db).split(), cwd=self.blastDataDir, env=os.environ, bufsize=2**16, stdout=PIPE, close_fds=True).stdout
        inp = readFastaRecords(pipe,readSeq=False)
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 500
        sql = """
        insert into seq
        (gi, taxid, src_db, project, seq_len, acc, kind, seq_hdr)
        values
        (%s,%s,%s,%s,%s,%s,%s,%s)
        """
        inserter = BulkInserter(cursor=curs,sql=sql,bufLen=bufLen)
        for rec in inp:
            (gifld,gi,accfld,acc,txt) = rec.title.split('|',4)
            assert gifld == 'gi'
            gi = int(gi)
            try:
                taxid = self.gi2taxa[gi]
            except KeyError:
                print "Taxid not found for gi: "+gi
                taxid = 0
            #Accession number formats are described here:
            #http://www.ncbi.nlm.nih.gov/Sequin/acc.html
            #WGS id might be 'AAAA00000000' or 'NZ_AAAA00000000', both cases seen in wgs db
            kind = ''
            acc_sfx = acc
            if acc_sfx[2] == '_':
                kind = acc_sfx[:2]
                acc_sfx = acc_sfx[3:]
            project = ''
            if acc_sfx[:4].isalpha() and acc_sfx[5].isdigit():
                project = acc_sfx[:4]
            values = (gi,taxid,db.id,project,rec.seqLen(),acc,kind,rec.title[:self.fastaHdrSqlLen])
            #values = [ str(x) for x in values ]
            inserter(values)            
            if iRec % 10000 == 0:
                print rec.title, taxid, iRec
            iRec += 1
        inserter.flush()
        curs.close()

    def loadRefseqAcc(self):
        self.dbSql.dropTable('refseq_acc')
        self.dbSql.ddl("""\
        create table refseq_acc
        (prefix varchar(2), acc varchar(40), molecule varchar(10), method varchar(15), descr varchar(400), index prefix(prefix))
        engine myisam
        """)
        self.dbSql.executemany("""\
        insert into refseq_acc
        (prefix,acc,molecule,method,descr)
        values
        (%s,%s,%s,%s,%s)
        """,refseqAccFormat)
        
    def loadTaxCategories(self):
        self.dbSql.dropTable('taxa_cat')
        self.dbSql.ddl("""\
        create table taxa_cat
        (cat char(1), taxid_ancestor integer, taxid integer)
        engine myisam
        """)
        inp = open(self.taxaCatFile,'r')
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 5000
        sql = """\
        insert into taxa_cat
        (cat,taxid_ancestor,taxid)
        values
        (%s,%s,%s)
        """
        inserter = BulkInserter(cursor=curs,sql=sql,bufLen=bufLen)
        for rec in inp:
            inserter(rec.split())
        inserter.flush()
        curs.close()
        self.dbSql.ddl("""\
        alter table taxa_cat  
        ADD PRIMARY KEY taxid(taxid), 
        add index cat(cat), 
        add index taxid_ancestor(taxid_ancestor)
        """)
        inp.close()

    def loadTaxNodes(self):
        self.dbSql.dropTable('taxa_node')
        self.dbSql.ddl("""\
        create table taxa_node
        (
        taxid integer,
        partaxid integer,
        rank varchar(20),
        embl_code  char(2),
        divid  integer,
        inh_div  bool,
        gcode_id integer,
        inh_gc  bool,
        mgcode_id integer,
        inhmgc  bool,
        gbhidden bool,
        hidsubtree bool,
        comments  varchar(40)
        )
        engine myisam
        """)
        inp = open(self.taxaNodesFile,'r')
        curs = self.dbSql.cursor()
        iRec = 0
        bufLen = 5000
        sql = """\
        insert into taxa_node
        (
        taxid,
        partaxid,
        rank,
        embl_code,
        divid,
        inh_div,
        gcode_id,
        inh_gc,
        mgcode_id,
        inhmgc,
        gbhidden,
        hidsubtree,
        comments
        ) 
        values
        (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
        """
        inserter = BulkInserter(cursor=curs,sql=sql,bufLen=bufLen)
        for rec in inp:
            inserter(rec.split('\t|\t'))
        inserter.flush()
        curs.close()
        #Get rid of spaces in 'rank' field to simplify export of data later on
        self.dbSql.execute("update taxa_node set rank = replace(rank,' ','_')").close()
        self.dbSql.ddl("""\
        alter table taxa_node
        ADD PRIMARY KEY taxid(taxid), 
        add index partaxid(partaxid), 
        add index divid(divid)
        """)
        inp.close()


    def loadTaxNodesMem(self):
        self.taxaTree = TaxaTree(ncbiDumpFile=self.taxaNodesFile)
        #self.taxaTree.write(sys.stdout)

        
    def loadGiTaxSql(self):
        self.dbSql.dropTable("gi_taxa")
        self.dbSql.execute(
        """
        create table gi_taxa
        (
        gi bigint,
        taxid integer
        )
        engine myisam
        """)
        #self.dbSql.executemany("insert into gi_tax (gi,taxid) values (%s,%s)",[(1,2),(3,4)])
        #print "Test done"
        #self.dbSql.execute("insert into gi_tax (gi,taxid) values ('1','2')")
        inp = file(self.taxaGiFile,'r')
        iRow = 0
        bufLen = 100000
        sql = "insert into gi_taxa (gi,taxid) values (%s,%s)"
        curs = self.dbSql.cursor()        
        inserter = BulkInserter(cursor=curs,sql=sql,bufLen=bufLen)        
        print "Start insert"
        for record in inp:
            inserter(record.split())
            iRow += 1
            if iRow % 100000 == 0:
                print record,iRow
        inserter.flush()
        curs.close()
        inp.close()
        print "Start index build"
        self.dbSql.ddl("""\
        alter table gi_taxa  
        ADD PRIMARY KEY gi(gi), 
        add index taxid(taxid)
        """)
        

    def loadGiTaxBdb(self,inpFile):
        import anydbm, whichdb
        print whichdb.whichdb('/export/atovtchi/taxa.db')
        gi2taxa = anydbm.open('/export/atovtchi/taxa.db', 'c')
        inp = file(inpFile,'r')
        iRow = 0
        buff = {}
        for line in inp:
            (gi,taxid) = line.split()
            buff[gi] = taxid
            if iRow % 100000 == 0:
                print gi, taxid, iRow
                gi2taxa.update(buff)
                buff = {}
            iRow += 1
        inp.close()
        taxaCnt = {}
        for (gi,taxid) in gi2taxa.iteritems():
            taxaCnt[taxid] = taxaCnt.get(taxid,0) + 1
        print sorted(((cnt,taxid) for (taxid,cnt) in taxaCnt.iteritems()))

    def loadGiTaxNumpy(self):
        gi2taxa = numpy.zeros(200*10**6,numpy.int)
        inp = file(self.taxaGiFile,'r')
        iRow = 0
        giMax = 0
        for rec in (line.split() for line in inp):
            gi2taxa[int(rec[0])] = int(rec[1])
            if iRow % 1000000 == 0:
                print rec, iRow
            giMax = max(giMax,int(rec[0]))
            iRow += 1
        inp.close()
        gi2taxa = numpy.resize(gi2taxa,giMax+1)
        out = file(self.taxaPickled,'w')
        dump(gi2taxa,out,-1)
        out.close()
        taxaCnt = {}
        for taxid in gi2taxa:
            taxaCnt[taxid] = taxaCnt.get(taxid,0) + 1
        print sorted(((cnt,taxid) for (taxid,cnt) in taxaCnt.iteritems()))

    def loadGiTaxPickled(self):
        inp = open(self.taxaPickled,'r')
        self.gi2taxa = load(inp)
        inp.close()


refseqAccFormat = \
(
('','','','','Undefined'),
('AC','AC_123456','Genomic','Mixed','Alternate complete genomic molecule. This prefix is used for records that are provided to reflect an alternate assembly or annotation. Primarily used for viral, prokaryotic records.'),
('AP','AP_123456','Protein','Mixed','Protein products; alternate protein record. This prefix is used for records that are provided to reflect an alternate assembly or annotation. The AP_ prefix was originally designated for bacterial proteins but this usage was changed.'),
('NC','NC_123456','Genomic','Mixed','Complete genomic molecules including genomes, chromosomes, organelles, plasmids.'),
('NG','NG_123456','Genomic','Mixed','Incomplete genomic region; supplied to support the NCBI genome annotation pipeline. Represents either non-transcribed pseudogenes, or larger regions representing a gene cluster that is difficult to annotate via automatic methods.'),
('NM','NM_123456 NM_123456789','mRNA','Mixed','Transcript products; mature messenger RNA (mRNA) transcripts.'),
('NP','NP_123456 NP_123456789','Protein','Mixed','Protein products; primarily full-length precursor products but may include some partial proteins and mature peptide products.'),
('NR','NR_123456','RNA','Mixed','Non-coding transcripts including structural RNAs, transcribed pseudogenes, and others.'),
('NT','NT_123456','Genomic','Automated','Intermediate genomic assemblies of BAC and/or Whole Genome Shotgun sequence data.'),
('NW','NW_123456 NW_123456789','Genomic','Automated','Intermediate genomic assemblies of BAC or Whole Genome Shotgun sequence data.'),
('NZ','NZ_ABCD12345678','Genomic','Automated','A collection of whole genome shotgun sequence data for a project. Accessions are not tracked between releases. The first four characters following the underscore (e.g. ''ABCD'') identifies a genome project.'),
('XM','XM_123456 XM_123456789','mRNA','Automated','Transcript products; model mRNA provided by a genome annotation process; sequence corresponds to the genomic contig.'),
('XP','XP_123456 XP_123456789','Protein','Automated','Protein products; model proteins provided by a genome annotation process; sequence corresponds to the genomic contig.'),
('XR','XR_123456','RNA','Automated','Transcript products; model non-coding transcripts provided by a genome annotation process; sequence corresponds to the genomic contig.'),
('YP','YP_123456 YP_123456789','Protein','Mixed','Protein products; no corresponding transcript record provided. Primarily used for bacterial, viral, and mitochondrial records.'),
('ZP','ZP_12345678','Protein','Automated','Protein products; annotated on NZ_ accessions (often via computational methods).')
)


    
