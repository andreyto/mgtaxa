### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

class CrisprGosApp(CrisprApp):
    """App-derived class for CRISPR array processing in GOS (Global Ocean Sampling) dataset"""

    BaseApp = CrisprApp
    
    def initWork(self,**kw):
        self.BaseApp.initWork(self,**kw)
        opt = self.opt
        self.mgtOutFile = "/usr/local/projects/GOSII/MGTAXA/v0.90/asm11_Mar06/genus/idlin.csv"
        self.aclameBlastFile = pjoin(opt.topWorkDir,"shannon/arr.blast.match_ACLAME.txt")
        self.aclameBlastOutFile = pjoin(opt.topWorkDir,"shannon/arr.blast.match_ACLAME.cov.txt")
    
    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "exportmgt":
            return self.exportMgtaxa()
        elif opt.mode == "exportaclame":
            return self.exportAclameBlast()
        elif opt.mode == "loadannot":
            return self.loadAnnot()
        elif opt.mode == "stats":
            return self.stats()
        else:
            return BaseApp.doWork(self,**kw)

    def exportMgtaxa(self):
        self.loadMgtaxa()

    def loadMgtaxa(self):
        """Load MGTAXA assignments into SQL tables.
        We parse the text output file, with each line like the one below:
        1118658405050   rank=genus : name=\"\"\"Candidatus Pelagibacter\"\"\" : taxid=198251 <<->>..."""
        dtype=[("id_seq","S%s"%self.maxSeqIdLen),
            ("taxid","i4"),
            ("rank","S20"),
            ("name","S80")]
        inp = openCompressed(self.mgtOutFile,'r')
        recs = []
        for line in inp:
            parts = line.split("taxid=",1)
            if len(parts) == 2:
                id_r = parts[0].split(None,1)
                id = 'scf'+id_r[0]
                rank = id_r[1].split('rank=')[1].split(':')[0].strip()
                name = id_r[1].split('name=')[1].split('"""')[1].strip()
                taxid = parts[1].split(None,1)[0]
                recs.append((id,taxid,rank,name))
        inp.close()
        recs = n.asarray(recs,dtype=dtype)
        db = self.getDbSql()
        db.createTableFromArray("arr_mgt",recs,indices=dict(primary="id_seq"))
        self.delDbSql()
    
    def loadAnnot(self):
        self.loadPepToRead()
        self.loadAsmToRead()
        self.loadPepAnnot()
        self.cleanCrisprPepAnnot()

    def loadPepToRead(self):
        inpFiles = [ "/usr/local/projects/GOSII/syooseph/clustering/phaseI_read_pred/site_phaseI_metagene_plus_clustering.gz",
                "/usr/local/projects/GOSII/syooseph/protein_predictions_GOS_IO/info_predictions.gz" ]
        dtype=[("id_pep","S%s"%self.maxSeqIdLen),
            ("id_read","S%s"%self.maxSeqIdLen)]
        recs = n.zeros(15000000,dtype=dtype)
        recsAp = ArrayAppender(recs)
        for inpFile in inpFiles:
            inp = openCompressed(inpFile,'r')
            for parts in ( line.split() for line in inp ):
                rec = (parts[0].split("JCVI_PEP_",1)[1].strip(), parts[5].split("JCVI_READ_",1)[1].strip())
                #for i in xrange(len(nextRec)): nextRec[i] = rec[i]
                arr,inext = recsAp.nextItem()
                arr[inext] = rec
            inp.close()
        recs = recsAp.getData()
        db = self.getDbSql()
        db.createTableFromArray("pep_read",recs,indices=dict(primary="id_pep",names=["id_read"]))
        self.delDbSql()
    
    def loadAsmToRead(self):
        inpFiles = [ "/usr/local/projects/GOSII/syooseph/asm11_Mar06_2009/mapping/asm_read_to_scf_deg_singleton_mapping.gz" ]
        dtype=[("id_asm","S%s"%self.maxSeqIdLen),
            ("begin","i8"),
            ("end","i8"),
            ("forward",bool),
            ("id_read","S%s"%self.maxSeqIdLen)]
        recs = n.zeros(15000000,dtype=dtype)
        recsAp = ArrayAppender(recs)
        for inpFile in inpFiles:
            inp = openCompressed(inpFile,'r')
            for parts in ( line.split() for line in inp ):
                rec = (parts[0],
                        int(parts[3]),
                        int(parts[4]),
                        parts[6] == '+',
                        parts[8].split("JCVI_READ_",1)[1].strip())
                # numpy record does not support slices [:]
                #(i.e. recsAp.nextElem()[:] = rec[:]) will raise InvalidIndex)
                # this will work:
                #nextRec = recsAp.nextElem()
                #for i in xrange(len(nextRec)): nextRec[i] = rec[i]
                arr,inext = recsAp.nextItem()
                arr[inext] = rec
            inp.close()
        recs = recsAp.getData()
        db = self.getDbSql()
        db.createTableFromArray("asm_read",recs,indices=dict(primary="id_read",names=["id_asm"]))
        self.delDbSql()

    def loadPepAnnot(self):
        annPat = "/usr/local/projects/GOSII/ANNOTATION/GS*/camera_annotation_rules/*/camera_annotation_rules.default.stdout.gz"
        dtype=[("id_pep","S%s"%self.maxSeqIdLen),
            ("name","S80")]
        recs = n.zeros(1000000,dtype=dtype)
        recsAp = ArrayAppender(recs)
        for inpFile in glob.glob(annPat):
            if not "GSIOVIR" in inpFile:
                print inpFile
                inp = openCompressed(inpFile,'r')
                for line in inp:
                    parts = line.split('\t')
                    #print '****'.join(parts)
                    id_pep = parts[0].split("JCVI_PEP_",1)[1].strip()
                    assert parts[1] == "common_name"
                    name = parts[2].strip()
                    if "crispr" in name.lower():
                        arr,inext = recsAp.nextItem()
                        arr[inext] = (id_pep,name)
                inp.close()
        recs = recsAp.getData()
        db = self.getDbSql()
        db.createTableFromArray("pep_ann",recs,indices=dict(primary="id_pep",names=["name"]))
        self.delDbSql()

    def cleanCrisprPepAnnot(self):
        db = self.getDbSql()
        recs = db.selectAsArray("""
        select *
        from pep_ann
        """)
        for rec in recs:
            name = rec["name"].item()
            name = name.replace("crispr","CRISPR")
            name = name.split("||")[0].strip()
            rec["name"] = name
        db.createTableFromArray("pep_ann_cr",recs,indices=dict(primary="id_pep",names=["name"]))
        self.delDbSql()
