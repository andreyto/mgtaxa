### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

class CrisprMooreApp(CrisprApp):
    """App-derived class for CRISPR array processing in "Moore genomes" dataset"""
    
    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "loadmic":
            return self.loadMic()
        elif opt.mode == "annotmic":
            return self.annotMicCrispr()
        elif opt.mode == "stats":
            return self.stats()
        else:
            return BaseApp.doWork(self,**kw)

    def loadMic(self):
        """Load data from Mic genomes necessary to run our CRISPR pipeline on them"""
        self.loadMooreSql()
        self.loadClosedSql()
        self.mergeMicPathsSql()
        self.makeMicNucFasta()

    def loadMooreSql(self):
        """Load meta-info about Moore genomes into SQL tables"""
        namesToScafFile = "/usr/local/projects/GOSII/syooseph/MF150/PROTEINS/final_init_len_gc"
        namesFile = "/usr/local/projects/GOSII/syooseph/MF150/PROTEINS/information_genome_len_gc_prot_tax.xls"
        dtype=[("path","S%s"%self.maxPathLen),
            ("name","S%s"%self.maxOrgNameLen),
            ("id_seq","S%s"%self.maxSeqIdLen)]
        inp = open(namesFile,'r')
        orgNames = set( ( parts[0] for parts in ( line.split('\t') for line in inp ) ) )
        inp.close()
        inp = open(namesToScafFile,'r')
        recs = [ (parts[0],parts[1],parts[2]) \
                for parts in ( line.split('\t') for line in inp ) \
                if parts[1] in orgNames ]
        # Glaciecola is special. We will need to load ALL loci from the corresponding file, which we indicate with '*'
        recs.append(("/usr/local/projects/GOSII/syooseph/MF150/Glaciecola_init/extracted_Glaciecola.gbf.gz",
            "Glaciecola sp. HTCC2999",
            "*"))
        inp.close()
        recs = n.asarray(recs,dtype=dtype)
        db = self.getDbSql()
        ## So far id_seq was unique in the dataset, at it is convenient to rely on it
        ## being unique in later processing. We define it as a primary key to be sure.
        db.createTableFromArray("moo_path",recs,indices=dict(primary="id_seq",names=["name"]))
        self.delDbSql()

    def loadClosedSql(self):
        """Load meta-info about Genbank non-Moore closed marine genomes into SQL tables"""
        namesToDirFile = "/usr/local/projects/GOSII/syooseph/MF150/MF_vs_closedgenomes/KEGG/parse_full_associate_KO/map_marine_orgs_names_proc"
        gbDirRoot = "/usr/local/projects/GOSII/syooseph/MF150/ClosedGenomes/data"
        dtype=[("path","S%s"%self.maxPathLen),
            ("name","S%s"%self.maxOrgNameLen),
            ("id_seq","S%s"%self.maxSeqIdLen)]
        inp = open(namesToDirFile,'r')
        orgNames =  [ (parts[0],parts[1]) \
                for parts in ( line.split('\t') for line in inp ) \
                if "closed" in parts[2].lower() ]
        inp.close()
        recs = []
        for (orgName,dirBase) in orgNames:
            dirName = pjoin(gbDirRoot,dirBase)
            drecs = []
            for gbFile in glob.iglob(pjoin(dirName,"*.gbk.gz")):
                drecs.append((gbFile,orgName,"*"))
            assert len(drecs) > 0, "No Genbank files found for organism %s in directory %s" % (orgName,dirName)
            recs += drecs
        recs = n.asarray(recs,dtype=dtype)
        db = self.getDbSql()
        db.createTableFromArray("clo_path",recs,indices=dict(primary="path,id_seq",names=["name"]))
        self.delDbSql()

    def mergeMicPathsSql(self):
        """Merge tables *_path with files paths for different sources of data, such as Moore and other Genbank marine genomes.
        @post Results are in mic_path table"""

        db = self.getDbSql()
        db.createTableAs("mic_path",
        """
        (
        SELECT  *,'moo' as src
        FROM    moo_path
        )
        """)
        
        db.ddl("""
        insert into mic_path
        (
        SELECT  *,'clo' as src
        FROM    clo_path
        )
        """)
        self.delDbSql()

    def makeMicNucFasta(self):
        """Export microbial GB nucl sequence into FASTA, store seq IDs into SQL table along the way, and split files for Piler run"""
        opt = self.opt
        db = self.getDbSql()
        pathArr = db.selectAsArray("select path,name,id_seq,src from mic_path")
        fileSeqs = groupRecArray(pathArr,"path")
        dtype=[("path","S%s"%self.maxPathLen),
            ("name","S%s"%self.maxOrgNameLen),
            ("id_seq","S%s"%self.maxSeqIdLen),
            ("src","S3"),
            ("genel","S3"),
            ("seq_len","i8")]
        del pathArr
        self.delDbSql()
        makedir(opt.crArrSeqDir)
        fastaOne = pjoin(opt.crArrSeqDir,self.crScafRoot+'.fasta.gz')
        out = openCompressed(fastaOne,'w')
        recs = []
        for fileName,recsFile in sorted(fileSeqs.items()):
            idSeqs = set(recsFile["id_seq"])
            orgName = list(set(recsFile["name"]))
            assert len(orgName) == 1, \
                    "Did not expect to see more than one organism name associated with a single Genbank file: %s %s" % (fileName,orgName)
            orgName = orgName[0]
            src = list(set(recsFile["src"]))[0]
            print "Reading GB file: %s" % fileName
            for seq in self.iterGenBankSeqs(fileName):
                # seq.name is seq.id w/o the '.version':
                if seq.name in idSeqs or "*" in idSeqs:
                    print orgName, seq.id, seq.description
                    genel = "unk"
                    if re.search(r'\Wplasmid\W',seq.description.lower()):
                        genel = "pla"
                    else:
                        for feat in seq.features:
                            if feat.type == "source":
                                if "plasmid" in feat.qualifiers:
                                    genel = "pla"
                                break
                    recs.append((fileName,orgName,seq.id,src,genel,len(seq.seq)))
                    SeqIO.write([seq],out,"fasta")
                else:
                    pdb.set_trace()
        out.close()
        
        recs = n.asarray(recs,dtype=dtype)
        db = self.getDbSql()
        db.createTableFromArray("mic_seq",recs,indices=dict(primary="id_seq",names=["name","path","src"]))
        self.delDbSql()

        makedir(opt.crArrDir)
        from bin.splitFasta import splitFasta
        splitFasta(inpFasta=fastaOne,outName=pjoin(opt.crArrDir,self.crScafRoot+'.fasta'),chunkSize=500*10**6)
    
    def annotMicCrispr(self):
        self.loadMicProtHitsSql()
        self.loadMicProtAnnotCrisprSql()
        self.annotMicCrisprGb()
        self.plotCrispr()
    
    def loadMicProtHitsSql(self):
        """Load reliable protein HMM hits for Moore and other closed genomes into SQL tables"""
        hitsFiles = ( "/usr/local/projects/GOSII/syooseph/MF150/HMM/processed/trusted_matches_info_genome.gz",
                "/usr/local/projects/GOSII/syooseph/MF150/ClosedGenomes/HMM_rest_marine/processed/trusted_matches_info_genome.gz",
                "/usr/local/projects/GOSII/syooseph/MF150/ClosedGenomes/HMM_Nitrosopumilus/processed/trusted_matches_info_genome.gz" )
        dtype=[("id","i4"),
            ("name_org","S%s"%self.maxOrgNameLen),
            ("id_prot","S%s"%self.maxSeqIdLen),
            ("id_model","S%s"%self.maxSeqIdLen),
            ("descr","S%s"%self.maxProtHitDescrLen),
            ("name_prot","S%s"%self.maxProtHitNameLen),
            ("class_prot","S%s"%self.maxProtHitNameLen)]
        recs = []
        id_rec = 0
        for hitsFile in hitsFiles:
            inp = openCompressed(hitsFile,'r')
            for rec in ( line.split('\t') for line in inp ):
                name_org = rec[0]
                id_prot = rec[1]
                if "|" in id_prot:
                    #always expect this, otherwise raise an exception
                    id_prot = re.split(r'\|(?:emb|gb|dbj|bbm)\|',id_prot)[1].split("|")[0]
                id_model = rec[2]
                descr = rec[3]
                class_prot = ""
                if "CRISPR" in descr:
                    name_prot = re.split(r":(?:.*\W)*CRISPR\W+",descr)
                    assert len(name_prot) == 2
                    name_prot = name_prot[0].upper()
                    if name_prot.startswith("CAS_"):
                        name_prot = name_prot.split("CAS_")[1]
                    elif name_prot.startswith("CRISPR"):
                        name_prot = name_prot.split("CRISPR_")[1]
                    class_prot = "CRISPR"
                elif "Crispr" in descr:
                    name_prot = descr.split(":")
                    assert len(name_prot) == 2
                    name_prot = name_prot[0]
                    class_prot = "CRISPR"
                if class_prot == "CRISPR":
                    recs.append((id_rec,name_org,id_prot,id_model,descr,name_prot,class_prot))
                    id_rec += 1
            inp.close()
        recs = n.asarray(recs,dtype=dtype)
        db = self.getDbSql()
        db.createTableFromArray("mic_prot_hits",recs,indices=dict(primary="id",names=["name_org","id_prot","name_prot","class_prot"]))
        self.delDbSql()

    def iterGenBankSeqs(self,gbFile):
        def _print_feat(feat):
            print str(feat)
            print feat.qualifiers.items()
        # some GB files have this "complement(1356723..&gt;1357840)", so we clean it up in a IO filter:
        inGb = ioFilter(openCompressed(gbFile,'r'),code="lambda x: x.replace('&gt;','').replace('&lt;','')",mode="line")
        for rec in SeqIO.parse(inGb,"genbank"):
            # for records with empty ID, set it to the value of the 'name' field
            if len(rec.id.strip())==0:
                rec.id = rec.name
            yield rec
        inGb.close()
    
    def loadMicProtAnnotCrisprSql(self):
        """Create SQL table with the most definitive name for CAS genes and their coordinates"""
        opt = self.opt
        dtype=[("id","i4"),
            ("name_org","S%s"%self.maxOrgNameLen),
            ("id_seq","S%s"%self.maxSeqIdLen),
            ("id_prot","S%s"%self.maxSeqIdLen),
            ("id_model","S%s"%self.maxSeqIdLen),
            ("name_prot","S%s"%self.maxProtHitNameLen),
            ("class_prot","S%s"%self.maxProtHitNameLen),
            ("begin","i8"),
            ("end","i8")]
        db = self.getDbSql()
        nameSeqs = groupRecArray(db.selectAsArray("""
        select name,path,id_seq from mic_seq"""),
        "name")
        protHits = groupRecArray(db.selectAsArray("""
        select * from mic_prot_hits
        """),
        "name_org")
        self.delDbSql()
        prot_recs = []
        prot_rec_id = 1
        for name,nameRecs in sorted(nameSeqs.items()):
            print "Processing organism: %s" % name
            fileSeqs = groupRecArray(nameRecs,"path")
            try:
                nameProtHits = groupRecArray(protHits[name],"id_prot")
            except KeyError:
                nameProtHits = {}
            for gbInp,recs in sorted(fileSeqs.items()):
                print "Reading GB file: %s" % gbInp
                for rec in self.iterGenBankSeqs(gbInp):
                    for feat in rec.features:
                        if feat.type == "CDS":
                            id_model = None
                            if "protein_id" in feat.qualifiers:
                                assert len(feat.qualifiers["protein_id"]) == 1
                                id_prot = feat.qualifiers["protein_id"][0]
                                jcvi_pref = "jcvi:cds."
                                if id_prot.startswith(jcvi_pref):
                                    id_prot_hit = id_prot.split(jcvi_pref)[1]
                                else:
                                    id_prot_hit = id_prot
                                if id_prot_hit in nameProtHits:
                                    pHits = nameProtHits.pop(id_prot_hit)
                                    best_pHit = [ pHit for pHit in pHits if pHit["id_model"].startswith("TIG") ]
                                    if len(best_pHit) > 0:
                                        best_pHit = best_pHit[0]
                                    else:
                                        best_pHit = pHits[0]
                                    name_prot = best_pHit["name_prot"]
                                    id_model = best_pHit["id_model"]
                                    print "Found HMM hit"
                                else:
                                    prod = feat.qualifiers.get("product",[""])[0].upper()
                                    if "CRISPR" in prod:
                                        name_prot = re.search(r"CRISPR-ASSOCIATED\W+(?:\w+\W+)*PROTEIN\W+(\w+)",prod)
                                        if name_prot:
                                            name_prot = name_prot.group(1)
                                        else:
                                            name_prot = re.search(r"CRISPR-ASSOCIATED\W+(.+)(?:\W+PROTEIN)*",prod)
                                            if name_prot is not None:
                                                name_prot = name_prot.group(1)
                                            else:
                                                name_prot = ''
                                                print "Could not parse out protein name"
                                        id_model = "ANNOT"
                                        print "Found product annotation in: %s" % (prod,)
                            else:
                                print "CDS feature does not have a protein_id qualifier: %s" % (feat,)

                            if id_model is not None:
                                prot_rec = (prot_rec_id,name,rec.id,id_prot,id_model,name_prot,"CRISPR",
                                    feat.location.nofuzzy_start,feat.location.nofuzzy_end)
                                prot_recs.append(prot_rec)
                                print "Adding protein record:\n%s" % (prot_rec,)
                                prot_rec_id += 1

        prot_recs = n.asarray(prot_recs,dtype=dtype)
        db = self.getDbSql()
        db.createTableFromArray("mic_prot_annot",prot_recs,indices=dict(primary="id",names=["name_org","id_seq","id_prot",
            "name_prot","class_prot","begin","end"]))
        self.delDbSql()

    def annotMicCrisprGb(self):
        """Add CRISPR arrays as annotation to Genbank format files"""
        opt = self.opt
        db = self.getDbSql()
        nameSeqs = groupRecArray(db.selectAsArray("""
        select name,path,id_seq from mic_seq 
        where 
        id_seq in 
        (select id_seq from arr) or 
        id_seq in
        (select id_seq from mic_prot_annot)"""),
        "name")
        arrs = groupRecArray(self.loadCrArrayRangesSql(db=db),"id_seq")
        arr_elems = groupRecArray(self.loadCrArrayElemSql(db=db),"id_arr")
        protAnn = groupRecArray(db.selectAsArray("""
        select * from mic_prot_annot
        """),
        "id_seq")
        self.delDbSql()
        makedir(opt.crAnnotDir)
        for name,nameRecs in sorted(nameSeqs.items()):
            print "Processing organism: %s" % name
            gbOut = pjoin(opt.crAnnotDir,strToFileName(name,remove=".,")+".gb")
            out = openCompressed(gbOut,'w')
            fileSeqs = groupRecArray(nameRecs,"path")
            for gbInp,recs in sorted(fileSeqs.items()):
                print "Reading GB file: %s" % gbInp
                for rec in self.iterGenBankSeqs(gbInp):
                    # Min range that encloses all features of interest.
                    # We expect CAS genes and CRISPR arrays to be close to each other
                    # We invert begin and end initial values here, so that 
                    # we can use min and max below
                    beginSub = len(rec)
                    endSub = 0
                    idSeqGb = rec.id
                    print rec.id, rec.description
                    nRecFeat = 0
                    if idSeqGb in protAnn:
                        recProts = protAnn[idSeqGb]
                        beginSub = recProts["begin"].min()
                        endSub = recProts["end"].max()
                        recProts = indexRecArray(recProts,"id_prot")
                        new_features = []
                        for feat in rec.features:
                            if feat.type == "CDS":
                                quals = feat.qualifiers
                                if "protein_id" in quals:
                                    assert len(quals["protein_id"]) == 1
                                    id_prot = quals["protein_id"][0]
                                    if id_prot in recProts:
                                        pAnn = recProts.pop(id_prot)
                                        feat.qualifiers["product"] = ("%s %s" % (pAnn["class_prot"],pAnn["name_prot"]),)
                                        print "Set protein qualifier"
                                        new_features.append(feat)
                                        nRecFeat += 1
                            else:
                                new_features.append(feat)
                        rec.features = new_features
                        assert len(recProts) == 0,"Some proteins were not found in SeqRecords" 
                    if idSeqGb in arrs:
                        recArrs = arrs[idSeqGb]
                        beginSub = min(beginSub,recArrs["begin"].min())
                        endSub = max(endSub,recArrs["end"].max())
                        for arr in recArrs:
                            id_ft = "CRISPR_%s" % arr["id_arr"]
                            ft = SeqFeature(location=FeatureLocation(arr["begin"],arr["end"]),
                                    type="repeat_region",
                                    strand=0,
                                    qualifiers={"rpt_type":("CRISPR",),"note":(id_ft,)},
                                    id=id_ft)
                            quals = ft.qualifiers
                            quals["rpt_unit_range"] = [ "%i..%i" % (el["begin"]+1,el["end"]) \
                                    for el in arr_elems[arr["id_arr"]] \
                                    if el["is_rep"] ]
                            rec.features.append(ft)
                            nRecFeat += 1
                    if beginSub < endSub:
                        # LOCUS name must be 16 or less:
                        # ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
                        # BioPython Genbank writer raises an exception when it has
                        # to use rec.name 'scf_xxxxxxxxxxx' as LOCUS name. For uniformity,
                        # we remove 'scf_' prefix from every rec.name that we output.
                        if rec.name.startswith('scf_'):
                            rec.name = rec.name[4:]
                        recSub = rec[max(beginSub-2000,0):endSub+2000]
                        assert len(recSub.features) >= nRecFeat,"Lost some features during slicing operation"
                        SeqIO.write([recSub],out,"genbank")
            out.close()


    def plotCrispr(self):
        """Use Bio.Graphics.GenomeDiagram to generate pictures of CRISPR genes and arrays from pre-processed Genbank files"""
        opt = self.opt
        rmdir(opt.crGraphDir)
        makedir(opt.crGraphDir)
        inpGbs = glob.glob(pjoin(opt.crAnnotDir,"*.gb"))
        for inpGb in inpGbs:
            outGraphFileRoot = pjoin(opt.crGraphDir,stripSfx(os.path.basename(inpGb)))
            for rec in self.iterGenBankSeqs(inpGb):
                self.plotCrisprOneGff(seqRec=rec,outGraphFileRoot=outGraphFileRoot)

    def plotCrisprOneBioGff(self,seqRec,outGraphFileRoot):
        """Generate GFF files of CRISPR genes and arrays from one SeqRecord from a pre-processed Genbank file, using Bio sandbox module"""

        from BCBio.GFF import GFF3Writer
        out = open(outGraphFileRoot+".gff3","w")
        writer = GFF3Writer()
        writer.write([seqRec],out)
        out.close()
    
    def plotCrisprOneGff(self,seqRec,outGraphFileRoot):
        """Generate GFF files and graphics for CRISPR genes and arrays from one SeqRecord from a pre-processed Genbank file"""

        from MGT.GFF import GFF3Record, GFF3Header
        from MGT.GFFTools import GFF3Graphics
        outFileRt = outGraphFileRoot+(".%s"%(seqRec.id,))
        gffFile = outFileRt + ".gff3"
        out = open(gffFile,"w")
        out.write(str(GFF3Header()))
        orec = GFF3Record(seqid=seqRec.id)
        for feat in seqRec.features:
            if feat.type in ("CDS","repeat_region"):
                for o in orec.fromSeqFeat(feat):
                    out.write(str(o))
        out.close()
        grFile = outFileRt + ".pdf"
        gr = GFF3Graphics(outFormat="pdf",width=max(len(seqRec)/10,800))
        gr(gffFile,grFile)
    
    def plotCrisprOneBioGd(self,seqRec,outGraphFileRoot):
        """Use Bio.Graphics.GenomeDiagram to generate pictures of CRISPR genes and arrays from one SeqRecord from a pre-processed Genbank file"""

        from reportlab.lib import colors
        from Bio.Graphics import GenomeDiagram

        gd_diagram = GenomeDiagram.Diagram(seqRec.id)
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
        gd_feature_set = gd_track_for_features.new_set()

        for feat in seqRec.features:
            if feat.type == "CDS" :
                color = colors.blue
                feat_name = feat.qualifiers["product"][0]
                if feat_name.upper().startswith("CRISPR"):
                    label_size = 14
                else:
                    label_size = 8
                gd_feature_set.add_feature(feat, sigil="ARROW",
                                           color=color, label=True,
                                           label_size = label_size, label_angle=30,
                                           name=feat_name)
            elif feat.type == "repeat_region" :
                color = colors.red
                gd_feature_set.add_feature(feat,
                                           color=color, label=True,
                                           label_size = 14, label_angle=30,
                                           name=feat.qualifiers.get("note",("repeat",))[0])
                
        #gd_diagram.draw(format="linear", pagesize='A4', fragments=4,
        #                start=0, end=len(seqRec))
        gd_diagram.draw(format="linear", pagesize=(3000,400), fragments=1,
                        start=0, end=len(seqRec))
        gd_diagram.write(outGraphFileRoot+".pdf", "PDF")
        gd_diagram.write(outGraphFileRoot+".svg", "SVG")

    def blastArrToGff(self,hsps,out):
        from MGT.GFF import GFF3Record, GFF3Header
        for hsp in hsps:
            orec = GFF3Record(seqid=hsp["id_q"],
                    Target="%s %s %s" % (hsp["id_h"],hsp["begin_h"]+1,hsp["end_h"]),
                    source="blast",
                    type="match_part",
                    score=hsp["mism"],
                    strand="+" if hsp["forward"]>0 else "-",
                    start=hsp["begin_q"],
                    end=hsp["end_q"]
                    )
            orec.attribs["Name"] = orec.attribs["Target"]
            out.write(str(orec))

