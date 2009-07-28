"""Horizontal chromosome transfer test.

1. Scan RefSeq nucleotide GenBank format files, and output multi-FASTA file
    a. Construct taxid_gi ID for each entry, followed by '| the original defline', skip plasmids
    b. Output all protein sequences (extracted from annotation features)
2. Concatenate the output of (1) for all microbes and for one diatom (to be used as outgroup)
3. Run cvtree on (2), which will produce the distance matrix
4. Run Phylip 'neighbour' on (3), which will produce the tree
5. Analyze and visualize the tree, looking for remotely placed chromosomes from the same taxid
    a. Color all multi-chromosome entries in red - there are few, we should see easily if some are not together
    b. Calculate tree distance between all entries from the same organism
"""


from MGT.Taxa import *
from MGT.FastaIO import FastaReader
from MGT.App import *

from MGT.DirStore import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from glob import iglob

class CvTreeId(Struct):
    """Class that encapsulates FASTA IDs we generate for CVTree programs"""
    @classmethod
    def splitFastaDeflineCvTree(klass,header):
        line = header.strip()
        if line.startswith('>'):
            #MGT.FastaReader
            line = line[1:]
        #BIO.SeqIO.parse() strips '>' itself
        parts = line.split('|')
        acc = parts[0]
        taxid,acc_gen = parts[1].split('_',1)
        return klass(acc=acc,taxid=int(taxid),acc_gen=acc_gen)

    @classmethod
    def fromOrgId(klass,orgid):
        self = klass()
        self.taxid,self.acc_gen = orgid.strip().split('_',1)
        self.taxid = int(self.taxid)
        return self

    def orgId(self):
        return "%s_%s" % (self.taxid,self.acc_gen)

class HctApp(App):
    """App-derived class for Horizontal chromosome transfer test"""

    #batchDepModes = ("parscan",)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.orgTypes = ["outgroup","microbial"]
        opt.maxOrgNameLen = 10
   
    def getDbSql(self):
        """Allocate (if necessary) and return a connection to SQL server"""
        if not hasattr(self,"dbSql"):
            self.dbSql = DbSqlMy(db="hct")
        return self.dbSql

    def delDbSql(self):
        """Call this to free a connection to SQL server if it will not be needed for extended period of time"""
        if hasattr(self,"dbSql"):
            self.dbSql.close()
            del self.dbSql

    def init(self):
        self.cvTreeExe = "/home/atovtchi/work/distros/cvtree/cvtree/cvtree"
        self.cvTreeMatExe = "/home/atovtchi/work/distros/cvtree/cvtree/batch_dist.pl"
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        self.cvSeqDir = self.store.getFilePath("cv-seq")
        self.cvVecDir = self.store.getFilePath("cv-vec")
        self.tmpDir = self.store.getFilePath("tmp")
        makedir(self.tmpDir)
        self.cvNamesFile = self.store.getFilePath("cv.name")
        self.cvMatFile = self.store.getFilePath("cv.mat")
        self.phMatFile = self.store.getFilePath("cv.phy.mat")
        self.idToNameFile = self.store.getFilePath("cv.id2name")
        self.treeDynAnnotFile = self.store.getFilePath("cv.tlf")
        self.outGroupAcc = 'NC_003424.3' # Schizosaccharomyces pombe 972h- chromosome I
        self.outGroupTaxid = 4896
        fungiDir = self.store.getFilePath("fungi")
        outGroupAccNoVer = stripSfx(self.outGroupAcc,'.')
        self.outGroupFna = pjoin(fungiDir,outGroupAccNoVer+'.fna')
        self.outGroupFaa = pjoin(fungiDir,outGroupAccNoVer+'.faa')
    
    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "seq":
            return self.pullSeq()
        elif opt.mode == "cv":
            return self.runCvTree()
        elif opt.mode == "mat":
            return self.buildDistMatrix()
        elif opt.mode == "phyinp":
            return self.renameMatrixTaxa()
        elif opt.mode == "cvmeta":
            return self.buildCvMetaData()
        elif opt.mode == "tree":
            return self.buildTree()

    def pullSeq(self):
        """Convert original RefSeq Fasta and GenBank files into SQL tables and protein Fasta files"""
        pass
        #self.loadGenomicIds()
        #self.pullProtSeq()
        #self.loadProtIds()
        #self.sqlPostPullSeq()
        #self.subSampleProtIds(maxProtSeqLen=100000)
        #protRecs = self.loadSelProtIdsMem()
        #self.makeCvTreeSeqInput(accProts=protRecs["acc"])
        #self.makeSeqGenNames()
        #self.exportTreeDynAnnot()
        #self.selectChromosomes(minProtSeqLen=100000)
        self.selectAnyGenElement(minProtSeqLen=100000)
        self.exportCvTreeNamesFile()

    def loadGenomicIds(self):
        db = self.getDbSql()
        self.createTableGenomicSeq()
        idGenSeq = IntIdGenerator()
        inserterSeq = db.makeBulkInserterFile(table='seq_gen',bufLen=50000,workDir=self.tmpDir)
        for orgType in self.opt.orgTypes:
            self.loadGenomicIdsOrgType(db=db,idGenSeq=idGenSeq,inserterSeq=inserterSeq,orgType=orgType)
        inserterSeq.flush()
        db.createIndices(table="seq_gen",
            names=["gi","taxid","acc","genel","orgtype"],
            primary="id")
        db.ddl("analyze table seq_gen",ifDialect="mysql")
        self.delDbSql()
    
    def loadGenomicIdsOrgType(self,db,idGenSeq,inserterSeq,orgType):
        """Load Fasta deflines for genomic sequences, index them by accession, and drop non-NC_ and all plasmids.
        @todo It might be more robust to parse the GenBank file, e.g.
        FEATURES             Location/Qualifiers
        source          1..208369
            /organism="Bacillus cereus ATCC 10987"
            /mol_type="genomic DNA"
            /strain="ATCC 10987"
            /db_xref="ATCC:10987"
            /db_xref="taxon:222523"
            /plasmid="pBc10987"
                                                                                                                           gene            join(207497..208369,1..687)

        We would have to fix the gap(unk100) bug first, and also check how the "extrachromosomal" is labeled
        in GB file.
        """
        if orgType == "outgroup":
            inFasta = self.outGroupFna
        else:
            inFasta = pjoin(options.refSeqDataDir,"%s.genomic.fna.gz" % orgType)
        inp = FastaReader(inFasta)
        for rec in inp.records():
            line = rec.header()[1:]
            parts = line.strip().split('|',4)
            assert parts[0] == 'gi'
            gi = int(parts[1])
            assert parts[2] == 'ref'
            acc = parts[3].strip()
            # we can do re.search(r'\bplasmid\b',) instead, but this is safer 
            # (there is a record called 'megaplasmid'):
            if acc[:3] == 'NC_':
                hdr = (' '.join(parts[4:])).strip()
                hlow = hdr.lower()
                ## genel values must differ by the first letter - this is used
                ## by the name generation method later
                if 'plasmid' in hlow:
                    genel = "pla"
                elif 'extrachromosomal' in hlow or 'extra-chromosomal' in hlow:
                    genel = "ext"
                elif 'transposon' in hlow:
                    genel = "tra"
                else:
                    genel = "chr"
                idSeq = idGenSeq()
                seqLen = rec.seqLen()
                values = (idSeq,gi,0,seqLen,acc,genel,orgType[:4],hdr)
                #values = [ str(x) for x in values ]
                inserterSeq(values)

        inp.close()

    def createTableGenomicSeq(self):
        db = self.getDbSql()
        db.ddl("""
        create table seq_gen
        (
        id integer,
        gi bigint,
        taxid integer,
        seq_len bigint,
        acc char(18),
        genel char(3),
        orgtype char(4),
        hdr varchar(80)
        )
        """,
        dropList=["table seq_gen"])

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def iterGenBankProts(self,gbFile):
        def _print_feat(feat):
            print str(feat)
            print feat.qualifiers.items()
        inGb = openCompressed(gbFile,'r')
        for rec in SeqIO.parse(inGb,"genbank"):
            #if iFoundRec > 20:
            #    break
            for feat in rec.features:
                quals = feat.qualifiers
                if feat.type == "source":
                    #_print_feat(feat)
                    taxid = int(dict([x.split(':') for x in quals['db_xref']])['taxon'])
                elif feat.type == "CDS":
                    #_print_feat(feat)
                    coded_by = quals["coded_by"][0]
                    #In [7]: re.findall(r"[A-Z]{2,}_\d+\.\d+",'complement(NC_006300.1:1934457..1934807)')
                    #Out[7]: ['NC_006300.1']
                    acc = re.findall(r"[A-Z]{2,}_\d+\.\d+",coded_by)[0]
            yield Struct(seq=rec.seq,accNa=acc,taxid=taxid,acc=rec.id)
        inGb.close()

    def pullProtSeq(self):
        for orgType in self.opt.orgTypes:
            if orgType == "outgroup":
                self.pullOutGroupProts()
            else:
                idsGen = self.loadGenomicIdsMem()
                self.pullGenBankProts(orgType=orgType,idsGen=idsGen)


    def pullGenBankProts(self,orgType,idsGen):
        """Scan a set of GenBank protein files and output a multi-FASTA file with each protein as a separate record and chromosome encoded in defline.
        The defline will look like:
        >YP_089573|221988_NC_006300.1
        as in ><protein_acc>|<taxid>_<chrmosome_acc>
        """
        gbFile    = pjoin(options.refSeqDataDir,"%s.protein.gpff.gz" % orgType)
        prFile = '%s.protein.cat.gz' % orgType
        outProt = openCompressed(prFile,'w')
        iRec = 0
        for rec in self.iterGenBankProts(gbFile):
            if rec.accNa in idsGen:
                seq = SeqRecord(rec.seq,
                        id = "%s|%s_%s" % (rec.acc,rec.taxid,rec.accNa),
                        description = '')
                SeqIO.write([seq],outProt,"fasta")
                iRec += 1
                if iRec % 100 == 0:
                    print "taxid = %s\t %s" % (rec.taxid,idsGen[rec.accNa])
        outProt.close()
            
    def pullOutGroupProts(self):
        """Convert RefSeq per-chromosome .faa file into .faa file with our defline format.
        The defline will look like:
        >YP_089573|221988_NC_006300.1
        as in ><protein_acc>|<taxid>_<chrmosome_acc>
        """
        prFile = 'outgroup.protein.cat.gz'
        outProt = openCompressed(prFile,'w')
        inProt = openCompressed(self.outGroupFaa,'r')
        for rec in SeqIO.parse(inProt,"fasta"):
            parts = rec.id.strip().split('|',4)
            assert parts[0] == 'gi'
            gi = int(parts[1])
            assert parts[2] == 'ref'
            acc = parts[3].strip()
            rec.id = "%s|%s_%s" % (acc,self.outGroupTaxid,self.outGroupAcc)
            rec.description = ''
            SeqIO.write([rec],outProt,'fasta')
        outProt.close()
        inProt.close()
            
    def pullGenBankProtsCatOld(self):
        """Scan a set of GenBank protein files and output a multi-FASTA file with all proteins in one chromosome concatenated.
        @todo We cannot insert spacers between proteins, because CvTree will not care. So we get chimeric k-mers."""
        def _out_prot():
            """Write the concat protein sequence"""
            if accNaLast in idsGen:
                recOut = SeqRecord(Seq(''.join(seq), rec.seq.alphabet),
                        id="%s_%s" % (taxidLast,accNaLast), 
                        description=idsGen[accNaLast])
                SeqIO.write([recOut],outProt,"fasta")
        #inGBs    = [ pjoin(options.refSeqDataDir,"microbial.genomic.gbff.gz") ]
        for orgType in self.opt.orgTypes:
            gbFile    = pjoin(options.refSeqDataDir,"%s.protein.gpff.gz" % orgType)
            idsGen = self.loadGenomicIds(orgType)
            prFile = '%s.protein.cat.gz' % orgType
            outProt = openCompressed(prFile,'w')
            iRec = 0
            iHostRec = 0
            iFoundRec = 0
            accNaLast = None
            taxidLast = None
            seq = []
            for rec in iterGenBankProts(gbFile):
                if accNaLast is None:
                    accNaLast = rec.accNa
                    taxidLast = rec.taxid
                if rec.accNa != accNaLast:
                    _out_prot()
                    seq = []
                    accNaLast = rec.accNa
                    taxidLast = rec.taxid
                seq.append(rec.seq.tostring())
            _out_prot()
            outProt.close()
            

    def loadProtIds(self):
        db = self.getDbSql()
        self.createTableProtSeq()
        idGenSeq = IntIdGenerator()
        inserterSeq = db.makeBulkInserterFile(table='seq_pro',bufLen=100000,workDir=self.tmpDir)
        for orgType in self.opt.orgTypes:
            self.loadProtIdsOrgType(db=db,idGenSeq=idGenSeq,inserterSeq=inserterSeq,orgType=orgType)
        inserterSeq.flush()
        db.createIndices(table="seq_pro",
            names=["taxid","acc","acc_gen","orgtype"],
            primary="id")
        db.ddl("analyze table seq_pro",ifDialect="mysql")
        self.delDbSql()
    
    def createTableProtSeq(self):
        db = self.getDbSql()
        db.ddl("""
        create table seq_pro
        (
        id integer,
        taxid integer,
        seq_len bigint,
        acc char(18),
        acc_gen char(18),
        orgtype char(4)
        )
        """,
        dropList=["table seq_pro"])
    
    def loadProtIdsOrgType(self,db,idGenSeq,inserterSeq,orgType):
        """Load Fasta deflines for protein sequences generated by pullSeq, parse and insert them into SQL table.
        """
        inFasta = self.store.getFilePath('%s.protein.cat.gz' % orgType)
        inp = FastaReader(inFasta)
        for rec in inp.records():
            hdr = CvTreeId.splitFastaDeflineCvTree(rec.header())
            idSeq = idGenSeq()
            seqLen = rec.seqLen()
            values = (idSeq,hdr.taxid,seqLen,hdr.acc,hdr.acc_gen,orgType[:4])
            inserterSeq(values)
        inp.close()

    def sqlPostPullSeq(self):
        """Create various derived tables after pulling sequences"""
        db = self.getDbSql()
        db.createTableAs("d_seq_pro_len","""
        (SELECT a.taxid,a.acc_gen,sum(a.seq_len) as seq_len
                FROM    seq_pro a
                GROUP BY taxid,acc_gen
        )
        """)
        db.createIndices(names=["taxid"],primary="acc_gen",table="d_seq_pro_len")
        ## Copy taxid from seq_pro to seq_gen
        db.ddl("""
        update seq_gen a, d_seq_pro_len b set a.taxid = b.taxid
        where a.acc = b.acc_gen
        """)
        ## Some seq_gen.taxid can stay zero after the above - for small plasmids with no annotated proteins
        db.createTableAs("d_seq_gen","""
        (SELECT a.*, b.seq_len as seq_pro_len
                FROM    seq_gen a, d_seq_pro_len b
                WHERE a.acc = b.acc_gen
        )
        """)
        db.createIndices(names=["taxid","acc","genel","gi","orgtype"],primary="id",table="d_seq_gen",
                attrib={"acc":{"unique":True},"gi":{"unique":True}})
        self.delDbSql()

    def makeSeqGenNames(self):
        """For each genetic element, we generate an artificial mnemonic name.
        The purpose of generating a new name is two-fold:
        1. Phylip cannot handle names longer than 10 symbols
        2. For each organism with multiple genetic elements, our name will reflect the
        nature of that element (chromosome, plasmid,...), as well as number elements of
        the same nature in the order of reducing size.
        Examples:
        If the organism is assigned id 123 and it has only one chromosome,
        that chromosome will be called t0123.
        If that organism has two chromosomes, they will be called t0123C1, and t0123C2,
        with C1 being the larger one.
        The disparity between one- and multi-chromosome naming schemes is designed to make
        multi-chromosomal organisms to visibly stand out.
        For plasmids, names will be like t0123P1.
        For all other genetic elements - like t0123?1 where ? is the first letter of 
        genel field of seq_gen table.
        The names are saved into an SQL table name_gen.
        """
        db = self.getDbSql()
        seqs = db.selectAsArray("""
        select taxid,acc,seq_len,genel
        from seq_gen
        where taxid <> 0
        order by taxid,genel,seq_len desc,acc
        """)
        seqs = groupRecArray(seqs,"taxid")
        names = []
        maxName = self.opt.maxOrgNameLen
        for (iOrg,taxid) in izipCount(sorted(seqs.keys())):
            orgName = "t%04i" % (iOrg + 1)
            assert len(orgName) <= maxName
            seqsOrg = groupRecArray(seqs[taxid],"genel")
            for (genel,seqGenel) in seqsOrg.items():
                if genel == "chr":
                    if len(seqGenel) == 1:
                        names.append((orgName,0,seqGenel[0]['acc'],orgName))
                        continue
                for (iEl,el) in izipCount(seqGenel):
                    name = "%s%s%s" % (orgName,genel[0].upper(),iEl+1)
                    assert len(name) <= maxName
                    names.append((orgName,1 if genel == "chr" else 0,el['acc'],name))
        accDtype = fieldDtypeRecArray(seqs.popitem()[1],'acc')
        names = n.asarray(names,dtype=[('orgname','S%s'%maxName),
            ('multichr',bool),
            ('acc',accDtype),
            ('name','S%s'%maxName)])
        db.createTableFromArray(name="d_seq_gen_name",arr=names)
        db.createIndices(names=["orgname","name"],primary="acc",table="d_seq_gen_name",
                attrib={"name":{"unique":True}})
        self.delDbSql() 

    def selectChromosomes(self,minProtSeqLen):
        """Select chromosomes subset for the distance matrix calculation"""
        db = self.getDbSql()
        db.createTableAs("d_seq_gen_sel","""
        (SELECT *
                FROM    d_seq_gen
                WHERE genel = "chr" and
                seq_pro_len >= %s and
                taxid <> 0

        )
        """ % (minProtSeqLen,))
        db.createIndices(names=["taxid","acc","genel","gi","orgtype"],primary="id",table="d_seq_gen_sel",
                attrib={"acc":{"unique":True},"gi":{"unique":True}})
        self.delDbSql()

    def selectAnyGenElement(self,minProtSeqLen):
        """Select genomic elements of any type longer than cutoff value for the distance matrix calculation"""
        db = self.getDbSql()
        db.createTableAs("d_seq_gen_sel","""
        (SELECT *
                FROM    d_seq_gen
                WHERE seq_pro_len >= %s and
                taxid <> 0

        )
        """ % (minProtSeqLen,))
        db.createIndices(names=["taxid","acc","genel","gi","orgtype"],primary="id",table="d_seq_gen_sel",
                attrib={"acc":{"unique":True},"gi":{"unique":True}})
        self.delDbSql()
    
    def exportCvTreeNamesFile(self):
        db = self.getDbSql()
        seqs = db.selectAsArray("""
        select taxid,acc
        from d_seq_gen_sel
        order by orgtype = 'outg' desc,taxid,acc
        """)
        self.delDbSql()
        names = [ CvTreeId(acc_gen=seq["acc"],taxid=seq["taxid"]).orgId() \
                for seq in seqs ]
        out = open(self.cvNamesFile,'w')
        out.write("%s\n" % len(names))
        out.write("\n".join(names)+"\n")
        out.close()

    def runCvTree(self):
        """Run cvtree executable that generates the composition vectors.
        Warning: this clears existing content of cvVecDir"""
        rmdir(self.cvVecDir)
        makedir(self.cvVecDir)
        run("%s -d %s -c %s -k 5 -t aa -S -i %s" \
                % (self.cvTreeExe,
                    self.cvSeqDir,
                    self.cvVecDir,
                    self.cvNamesFile))

    def runCvTreeMat(self):
        """Run cvtree script that builds the distance matrix from precomputed composition vectors.
        """
        rmrf(self.cvMatFile)
        run("%s 2 %s %s %s" \
                % (self.cvTreeMatExe,
                    self.cvNamesFile,
                    self.cvVecDir,
                    self.cvMatFile),
                shell=True)

    def buildDistMatrix(self):
        self.runCvTreeMat()
        self.cvMatToPhylipMat()

    def subSampleProtIds(self,maxProtSeqLen):
        """Select a subset of proteins from each genomic sequence constrained by a sum of protein lengths.
        The use case: we want to have each chromosome represented by an equal length sample, to make
        sure that the phylogenetic tree we build is not influenced by length differences."""
        db = self.getDbSql()
        prots = db.selectAsArray("""
        select acc_gen,acc,seq_len 
        from seq_pro
        order by acc_gen,acc
        """)
        genProts = groupRecArray(prots,"acc_gen")
        selRecs = []
        for acc_gen in genProts.keys():
            recs = nrnd.permutation(genProts[acc_gen])
            indMax = recs["seq_len"].cumsum().searchsorted(maxProtSeqLen,side="right")
            selRecs.append(recs[:indMax])
        selRecs = n.concatenate(selRecs)
        db.createTableFromArray(name="d_seq_pro_sel_1",arr=selRecs)
        self.delDbSql()

    def loadSelProtIdsMem(self):
        """Return as Numpy array table created by subSampleProtIds()"""
        db = self.getDbSql()
        ret = db.selectAsArray("""
        select * 
        from d_seq_pro_sel_1
        order by acc_gen,acc
        """)
        self.delDbSql()
        return ret


    def makeCvTreeSeqInput(self,accProts):
        """Create input FASTA files for CVTree by selecting only sequences with Acc from accProts.
        Warning: this will erase existing content of the target directory."""
        rmdir(self.cvSeqDir)
        makedir(self.cvSeqDir)
        accProts = set(accProts)
        acc_gen_null = 'XXXXXXXXXXXXX'
        acc_gen_last = acc_gen_null
        for orgType in self.opt.orgTypes:
            print "Processing sequences of orgType '%s'" % orgType
            #idsGen = self.filterGenomicIds(self.loadGenomicIds(orgType))
            prFile = '%s.protein.cat.gz' % orgType
            inProt = openCompressed(prFile,'r')
            for rec in SeqIO.parse(inProt,"fasta"):
                hdr = CvTreeId.splitFastaDeflineCvTree(rec.id)
                if hdr.acc in accProts:
                    #print rec.id
                    if hdr.acc_gen != acc_gen_last:
                        if acc_gen_last != acc_gen_null:
                            out.close()
                        acc_gen_last = hdr.acc_gen
                        #CvTree expects '.faa' extension for proteins
                        outName = pjoin(self.cvSeqDir,hdr.orgId()+'.faa')
                        assert not os.path.exists(outName),"We expect proteins grouped by genetic element in the input FASTA file"
                        out = open(outName,'w')
                    rec.description = ''
                    SeqIO.write([rec],out,'fasta')
        if acc_gen_last != acc_gen_null:
            out.close()


    def listAllCvTreeNames(self):
        return [ stripSfx(os.path.basename(f),'.faa') for f in iglob(pjoin(self.cvSeqDir,'*.faa')) ]

    def loadCvTreeNamesFile(self):
        """Load names from a name file"""
        inp = open(self.cvNamesFile,'r')
        l = [ y for y in ( x.strip() for x in inp ) if len(y) > 0 ]
        inp.close()
        return l
    
    def cvMatToPhylipMat(self):
        """Rename sequence names in the distance matrix into Phylip compatible names"""
        accToName = self.loadAccToNameMem()
        inpMat = open(self.cvMatFile,'r')
        outMat = open(self.phMatFile,'w')
        line = inpMat.next()
        outMat.write(line)
        for line in inpMat:
            cvName,rest = line.split(None,1)
            acc = CvTreeId.fromOrgId(cvName).acc_gen
            phName = accToName[acc]
            assert len(phName) <= 10
            outMat.write("%s%s" % (phName.ljust(10),rest))
        outMat.close()
        inpMat.close()

    def loadAccToNameMem(self):
        """Load dict(acc->name) from SQL table"""
        db = self.getDbSql()
        ret = db.selectAsNx1Dict("""
        select acc,name 
        from d_seq_gen_name
        """)
        self.delDbSql()
        return ret

    def loadSeqGenAnnotMem(self):
        db = self.getDbSql()
        ret = db.selectAsArray("""
        select b.name as id, 
        a.taxid,a.acc,
        a.genel,a.orgtype,
        a.seq_len,a.seq_pro_len,
        b.multichr
        from d_seq_gen a, d_seq_gen_name b
        where a.acc = b.acc
        order by b.name
        """)
        self.delDbSql()
        return ret

    
    def exportTreeDynAnnot(self):
        """Create a file in TreeDyn annotation format.
        Example of such file format (taken from TreeDyn web site):
        BUD2  Subcellular_loc { Bud_neck Cytoskeletal } Cellular_Role { Cell_polarity } FuncCat { GTPase_activating_protein }
        BUD3  Subcellular_loc { Bud_neck } Cellular_Role { Cell_polarity } FuncCat { Unknown }
        """
        def _quote_taxa_name(tname):
            return tname.replace(' ','_')[:20]
        def _format_non_id_fields(rec):
            return ' '.join(["%s { %s }" % (fld_name,rec[fld_name]) 
                    for fld_name in sorted(rec.dtype.fields.keys())
                    if fld_name != 'id'])

        recs = self.loadSeqGenAnnotMem()
        taxaTree = self.getTaxaTree()
        taxaLevels = TaxaLevels(taxaTree)
        levelNames = taxaLevels.getLevelNames()
        out = open(self.treeDynAnnotFile,"w")
        for rec in recs:
            taxid = rec["taxid"]
            node = taxaTree.getNode(taxid)
            lin = taxaLevels.lineageFixedList(node)
            lin_s = ' '.join(["%s { %s }" % (lev_n,
                _quote_taxa_name(taxaTree.getNode(lev_v).name) if lev_v is not None else None) 
                for (lev_n,lev_v) in it.izip(levelNames,lin)])
            line = "%s %s %s\n" % (rec["id"],_format_non_id_fields(rec),lin_s)
            out.write(line)
        out.close()


def run_Hct():
    opt = Struct()
    opt.runMode = "inproc" #"batchDep"
    modes = ["seq","cv","mat"] #"mat" "cv" "seq"
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = HctApp(opt=opt)
        jobs = app.run(depend=jobs)

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(HctApp)

