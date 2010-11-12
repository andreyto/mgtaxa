### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Class to process CRISPR arrays found by tools like PILERCR"""

from MGT.Taxa import *
from MGT.FastaIO import FastaReader
from MGT.BlastDb import BlastDb,runBlast
from MGT.App import *

from MGT.DirStore import *

from MGT import UUID

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio import pairwise2
from Bio.Blast import NCBIXML


def crisprLengthRatioOk(arr,spacerToRepeatRatio=(0.6,2.5)):
    spLen = n.asarray([len(rec['spacer']) for rec in arr[:-1]])
    repLen = n.asarray([len(rec['repeat']) for rec in arr])
    avgRepLen = repLen.mean()
    if not n.logical_and(spLen >= avgRepLen*spacerToRepeatRatio[0],spLen <= avgRepLen*spacerToRepeatRatio[1]).all():
        return False
    return True

def crisprHasSimilarSpacers(arr,nLookAhead=3,maxIdentity=0.6,maxSimilarRatio=0.2):
    aligner = pairwise2.align.globalxs
    # Copy of comments from Bio.pairwise2 code:
    # The alignment functions take some undocumented keyword parameters:
    # - penalize_extend_when_opening: boolean
    #   Whether to count an extension penalty when opening a gap.  If
    #   false, a gap of 1 is only penalize an "open" penalty, otherwise it
    #   is penalized "open+extend".
    # - penalize_end_gaps: boolean
    #   Whether to count the gaps at the ends of an alignment.  By
    #   default, they are counted for global alignments but not for local
    #   ones.
    # - gap_char: string
    #   Which character to use as a gap character in the alignment
    #   returned.  By default, uses '-'.
    # - force_generic: boolean
    #   Always use the generic, non-cached, dynamic programming function.
    #   For debugging.
    # - score_only: boolean
    #   Only get the best score, don't recover any alignments.  The return
    #   value of the function is the score.
    # - one_alignment_only: boolean
    #   Only recover one alignment.
    spacers = arr[:-1]
    found = 0
    for i in xrange(len(spacers)):
        for j in xrange(i+1,min(i+1+nLookAhead,len(spacers))):
            alment = aligner(\
                    spacers[i]['spacer'],
                    spacers[j]['spacer'],
                    -0.8,-0.2,penalize_end_gaps=1,one_alignment_only=0)
            align1, align2, score, begin, end = alment[0]
            #if spacers[i]['pos'] == 8627:
            #    pdb.set_trace()
            a1 = n.fromstring(align1,dtype='S1')
            a2 = n.fromstring(align2,dtype='S1')
            nident = (a1 == a2).sum()
            ident = float(nident)/len(a1)
            if ident >= maxIdentity:
                found+=1
                break
    return found > maxSimilarRatio * len(spacers)

def crisprSeqAliToRaw(seq_ali):
    return seq_ali.replace("-","")

def printAliSeqs(seqs,lineLen,out,seqNames=None,emptySymb=' '):
    maxLen = max((len(s) for s in seqs))
    if seqNames is None:
        seqNames = ['']*len(s)
    seqs = [ s.ljust(maxLen,emptySymb) for s in seqs ]
    maxNameLen = max((len(s) for s in seqNames))
    for x in range(0,maxLen,lineLen):
        for (s,sName) in it.izip(seqs,seqNames):
            out.write(sName.ljust(maxNameLen+1))
            out.write(s[x:x+lineLen])
            out.write("\n")
        out.write("\n")

class RangeSet(object):
    """Class that converts ranges (same as intervals) into sets in order to find intersections and coverage"""
    
    def __init__(self):
        self.s = set()

    def addRanges(self,begin,end):
        """Add several ranges represented by begin and end sequences"""
        s = self.s
        for (b,e) in it.izip(begin,end):
            s.update(xrange(b,e))

    def toSet(self):
        return self.s


class CrisprApp(App):
    """App-derived base class for CRISPR array processing"""

    batchDepModes = ("blastcr","pilercr","findcr")

    BaseApp = App

    maxSeqIdLen = UUID.maxIdLen
    maxSeqPartIdLen = UUID.maxIdLen
    maxArrElemLen = 200
    maxPathLen = 200
    maxOrgNameLen = 100
    maxProtHitDescrLen = 100
    maxProtHitNameLen = UUID.maxIdLen

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        globOpt = globals()["options"]
        opt = options
        opt.setIfUndef("sqlHost","mgtaxa-dev.jcvi.org")
        opt.setIfUndef("sqlDb","crispr")
        opt.minSpacerLen = 20
        opt.maxSpacerLen = 60
        opt.minSpacerNum = 3
        workDir = opt.setIfUndef("topWorkDir",globOpt.dataDir) #os.environ["GOS_WORK"]
        opt.crArrSeqDir = pjoin(workDir,"scaf-crispr")
        opt.crArrDir = pjoin(workDir,"scaf-piler")
        opt.gosBlastDbDir = pjoin(workDir,"reads-blast")
        opt.crBlastDir = pjoin(workDir,"crispr-blast")
        opt.crAnnotDir = pjoin(workDir,"crispr-annot")
        opt.crGraphDir = pjoin(workDir,"crispr-graph")
        opt.MEM = 6000
        ## Combine all viral tracks and all microbial tracks for BLAST array hits
        opt.virMicTracksUnion = True
        ## Only output tracks if viral hits are present for BLAST array hits
        opt.virArrHitsOnly = True
   
    def getDbSql(self):
        """Allocate (if necessary) and return a connection to SQL server"""
        if not hasattr(self,"dbSql"):
            self.dbSql = DbSqlMy(db=self.opt.sqlDb,host=self.opt.sqlHost)
        return self.dbSql

    def delDbSql(self):
        """Call this to free a connection to SQL server if it will not be needed for extended period of time"""
        if hasattr(self,"dbSql"):
            self.dbSql.close()
            del self.dbSql

    def initWork(self,**kw):
        self.BaseApp.initWork(self,**kw)
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        makedirs((opt.topWorkDir,opt.crArrSeqDir,opt.crArrDir,opt.crBlastDir,opt.crAnnotDir,opt.crGraphDir))
        self.tmpDir = pjoin(opt.topWorkDir,"tmp")
        makedir(self.tmpDir)
        storeBlast = SampStore.open(path=opt.crBlastDir)
        self.crArrSeqFile = storeBlast.getFilePath("arr.fasta")
        self.crArrBlastFileRoot = storeBlast.getFilePath("arr.blast")
        self.crArrBlastTxtFile = storeBlast.getFilePath("arr.blast.txt")
        self.crArrBlastPerMatchTxtFile = storeBlast.getFilePath("arr.blast.match.txt")
        self.crArrBlastMatchFile = storeBlast.getFilePath("arr.blast.match.fasta")
        self.initBlastJobs()
        self.crScafRoot = "inp"
   
    def initBlastJobs(self):
        opt = self.opt
        self.pandaBlastDbDir = "/usr/local/db/panda/NuclSequences"
        pandaBlastGroups = ["BacterialGroup", "ArchaeaGroup", "PlasmidGroup", "ViralGroupNucl", "Phage"]
        gosBlastGroups = [ "vir" ]
        pandas = [ Struct(dbPath=pjoin(self.pandaBlastDbDir,gr+".fasta"),
            outFile=self.crArrBlastFileRoot+".panda."+gr,
            inpFile=self.crArrSeqFile,
            group = gr,
            kind = "panda",
            order = ord + 10) for \
                (ord,gr) in izipCount(pandaBlastGroups) ]
        goses = [ Struct(dbPath=pjoin(opt.gosBlastDbDir,gr),
            outFile=self.crArrBlastFileRoot+".gos."+gr,
            inpFile=self.crArrSeqFile,
            group = gr,
            kind = "gos",
            order = ord + 1) for \
                (ord,gr) in izipCount(gosBlastGroups) ]
        self.blastDbs = pandas + goses
        for db in self.blastDbs:
            db.outFileBin = db.outFile+".pkl.gz"

    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "findcr":
            return self.findCr(**kw)
        elif opt.mode == "pilercr":
            return self.pilerCr(**kw)
        elif opt.mode == "pilercr-one":
            return self.pilerCrOne(**kw)
        elif opt.mode == "loadcr":
            return self.loadCrispr()
        elif opt.mode == "exportcr":
            return self.exportCr()
        elif opt.mode == "blastdb":
            return self.makeBlastDb()
        elif opt.mode == "blastcr":
            return self.blastCr(**kw)
        elif opt.mode == "blastcr-one":
            return self.blastCrOne(**kw)
        elif opt.mode == "loadbl":
            return self.loadBlast()
        elif opt.mode == "stats":
            return self.stats()
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def findCr(self,**kw):
        """Run full pipeline for CRISPR loci detection - from finder to SQL load and FASTA export"""
        opt = self.opt.copy()
        opt.mode = "pilercr"
        app = self.factory(opt=opt)
        jobs = app.run(**kw)
        opt = self.opt.copy()
        opt.mode = "loadcr"
        app = self.factory(opt=opt)
        jobs = app.run(depend=jobs)
        opt = self.opt.copy()
        opt.mode = "exportcr"
        app = self.factory(opt=opt)
        jobs = app.run(depend=jobs)
        return jobs


    def loadCrispr(self):
        """Import and filter original PILERCR output files and convert them into SQL tables"""
        arr = self.importPiler()
        self.pilerToSql(arr)

    def exportCr(self):
        self.exportCrisprArraySeq()

    def makeBlastDb(self):
        opt = self.opt
        curdir = os.getcwd()
        makedir(opt.gosBlastDbDir)
        os.chdir(opt.gosBlastDbDir)
        # formatdb take a single space separatated string for -i option
        # when you need to combine several FASTA files to a single DB,
        # and it is crazily sensitive to delimiteres, that is why find -printf is used below
        # if we need to search in subdirectories.
        #seqFiles=$(find GSIO* -name inp.fasta -printf " %p")
        #seqFiles=backsticks("""find . -name '*.fasta' -printf " %p" """,shell=True)
        from bin.gos import makeReads
        seqFiles = makeReads.listReadFiles()
        gosBlastDb = os.path.basename([ db.dbPath for db in self.blastDbs if db.kind == "gos" ][0])
        run("""formatdb -n %s -t %s -p F -i "%s" """ % (gosBlastDb,gosBlastDb,' '.join(seqFiles)),shell=True,debug=1)
        os.chdir(curdir)

    def blastCr(self,**kw):
        """BLAST CRISPR arrays against different DBs in parallel"""
        opt = self.opt
        #pandaDbs = ["Phage.fasta"]
        jobs = []
        for blDb in self.blastDbs:
            blOpt = copy(opt)
            blOpt.mode = "blastcr-one"
            blOpt.blDbPath = blDb.dbPath
            blOpt.blOutFile = blDb.outFile
            blOpt.blOutFileBin = blDb.outFileBin
            blOpt.blInpFile = blDb.inpFile
            blOpt.MEM = 6000
            blOpt.LENGTH = "medium"
            blOpt.ARCH = "lx*64"
            blApp = self.factory(opt=blOpt)
            print blApp.opt
            jobs.extend(blApp.run(**kw))
        return jobs

    def blastCrOne(self,**kw):
        """One job launched by blastCr()"""
        opt = self.opt
        runBlast(dbPath=opt.blDbPath,inpFile=opt.blInpFile,outFile=opt.blOutFile,paramStr="-p blastn -e 100 -W 7 -F F -m 7 -U F")
        # Now parse with relaxed cutoffs and save as numpy record array for faster loading later
        recs = self.parseArrayBlast(inFile=opt.blOutFile,minAlignLen=10,minBitScore=10,maxEvalue=100,maxMism=100,out=None,debug=False)
        dumpObj(recs,opt.blOutFileBin)

    def pilerCr(self,**kw):
        """Find CRISPR arrays with PilerCr within several partitions of the input sequence in parallel"""
        opt = self.opt
        inpFiles = self.listPilerIoFiles(input=True)
        jobs = []
        for inpFile in inpFiles:
            jOpt = copy(opt)
            jOpt.mode = "pilercr-one"
            jOpt.plInpPath = inpFile
            jOpt.MEM = 2000
            jOpt.LENGTH = "medium"
            jOpt.ARCH = "lx*64"
            jApp = self.factory(opt=jOpt)
            print jApp.opt
            jobs.extend(jApp.run(**kw))
        return jobs

    def pilerCrOne(self,**kw):
        """One job launched by pilerCrOne()"""
        opt = self.opt
        plOutPath = pjoin(opt.crArrDir,os.path.basename(opt.plInpPath)+".piler")
        run(("/home/atovtchi/work/distros/CRISPR/PILERCR/src/pilercr -minrepeat 22 -in %s -out %s -outtab %s.csv" % \
                (opt.plInpPath,plOutPath,plOutPath)).split())
    
    def blastToCoverage(self):
        """Calculate how much spacers vs repeats each BLAST match hits.
        @return record array with (id_arr,id_hit,coverage_spacers,coverage_repeats)"""
        hsps = []
        for blDb in self.blastDbs:
            if blDb.group in ("vir",):
                hspsOne = self.loadArrayBlastBin(blDb.outFileBin,minAlignLen=22,minBitScore=1,maxEvalue=1000,maxMism=2,out=None)
                hsps.append(hspsOne)
        hsps = n.concatenate(hsps)
        hsps = groupRecArray(hsps,"id_q")
        arr = self.loadCrisprSqlArr()
        arr = groupRecArray(arr,"id_arr")
        dtype=[("id_q","i4"),
            ("id_h","S%s"%self.maxSeqIdLen),
            ("cov_sp","i4"),
            ("cov_re","i4")]
        cov = []
        for (id_arr,arr_els) in sorted(arr.items()):
            # Array coords are relative to scaffolds, but we aligned to extracted array sequences,
            # so the HSP coords are relative to the start of the array. We modify the array coords
            # assuming that self.loadCrisprSqlArr() returns array elements in coord order.
            # @todo maybe we should instead convert to scaffold coords when we load Blast hits?
            beg = arr_els[0]["begin"]
            arr_els["begin"] -= beg
            arr_els["end"] -= beg
            if id_arr in hsps:
                # Make range sets separately for all spacers and all repeats in array
                ra_sp = RangeSet()
                arr_els_sp = arr_els[arr_els["is_rep"] == 0]
                ra_sp.addRanges(arr_els_sp["begin"],arr_els_sp["end"])
                ra_re = RangeSet()
                arr_els_re = arr_els[arr_els["is_rep"] != 0]
                ra_re.addRanges(arr_els_re["begin"],arr_els_re["end"])
                # Loop through all hsps that hit this array, grouped by match id
                arr_hsps = groupRecArray(hsps[id_arr],"id_h")
                for (id_h,hsps_h) in sorted(arr_hsps.items()):
                    # Make range set for hsps coords on the array
                    ra_h = RangeSet()
                    ra_h.addRanges(hsps_h["begin_q"],hsps_h["end_q"])
                    cov_sp = len(ra_sp.toSet()&ra_h.toSet())
                    cov_re = len(ra_re.toSet()&ra_h.toSet())
                    cov.append((id_arr,id_h,cov_sp,cov_re))
        return n.asarray(cov,dtype=dtype)
                    

    def exportAclameBlast(self):
        cov = self.blastToCoverage()
        cov = groupRecArray(cov,"id_h")
        inp = openCompressed(self.aclameBlastFile,'r')
        out = openCompressed(self.aclameBlastOutFile,'w')
        lines = inp.readlines()
        if len(lines) == 1:
            lines = lines[0].split('\r')
        inp.close()
        for line in lines:
            if line.startswith('#'):
                line = '#'+'\t'.join(["id_array","coverage_spacers","coverage_repeats"])+'\t'+line[1:]+'\n'
                out.write(line)
            else:
                parts = line.strip().split('\t')
                id_h = parts[1].strip()
                cov_h = cov[id_h]
                assert len(cov_h) > 0
                for c in cov_h:
                    parts_o = copy(parts)
                    parts_o = [ "%s" % x for x in ( c["id_q"],c["cov_sp"],c["cov_re"] ) ] + parts_o
                    out.write('\t'.join(parts_o)+'\n')
        out.close()
            

    def loadBlast(self):
        opt = self.opt
        hspsTracks = {}
        for blDb in self.blastDbs:
            #if blDb.group not in ("Phage",):
            #    continue
            hsps = self.loadArrayBlastBin(blDb.outFileBin,minAlignLen=22,minBitScore=1,maxEvalue=1000,maxMism=1,out=None)
            if opt.virMicTracksUnion:
                if blDb.group in ("ViralGroupNucl","Phage", "vir"):
                    trackName  = "01.vir"
                else:
                    trackName = "10.mic"
            else:
                trackName = ("%2i"%blDb.order)+'.'+blDb.kind+'.'+blDb.group
            if hspsTracks.has_key(trackName):
                hspsTracks[trackName] = n.concatenate([hspsTracks[trackName],hsps])
            else:
                hspsTracks[trackName] = hsps
        seqsQ = self.loadArraySeqFile()
        self.exportBlast(seqsQ=seqsQ,hspsTracks=hspsTracks)
        #self.exportBlastPerMatch(seqsQ=seqsQ,hspsTracks=hspsTracks)
        #self.exportBlastMatches(hsps=hsps)

        from MGT.GFF import GFF3Record, GFF3Header
        out = open("tmp.gff3","w")
        out.write(str(GFF3Header()))
        self.blastArrToGff(hspsTracks["01.vir"],out)
        out.close()

    def exportBlast(self,seqsQ,hspsTracks):
        opt = self.opt
        res = defdict(Struct)
        for (trackName,hsps) in sorted(hspsTracks.items()):
            hsps = groupRecArray(hsps,"id_q")
            for (id_q,hsps_q) in sorted(hsps.items()):
                seq_q = seqsQ[id_q]
                seq_a_t = n.zeros(len(seq_q),dtype='S1')
                seq_a_t[:] = '.'
                #out.write("\n%s\n" % seq_q)
                seq_a = seq_a_t.copy()
                for hsp in hsps_q:
                    seq_ali_h = n.fromstring(hsp["seq_ali_h"],dtype='S1')
                    try:
                        seq_a[hsp["begin_q"]:hsp["end_q"]] = seq_ali_h
                    except:
                        ali_max_len = fieldDtypeRecArray(hsps_q,"seq_ali_h").itemsize
                        ali_q_len = -hsp["begin_q"] + hsp["end_q"]
                        if ali_q_len > ali_max_len:
                            print "Alignment truncated due to limited field length"
                            seq_a[hsp["begin_q"]:hsp["begin_q"]+ali_max_len] = seq_ali_h
                        else:
                            print "Unexplained length mismatch in alignment, truncating: %s %s" % (ali_q_len,len(seq_ali_h))
                            sublen = min(len(seq_ali_h),len(seq_a)-hsp["begin_q"])
                            seq_a[hsp["begin_q"]:hsp["begin_q"]+sublen] = seq_ali_h[:sublen]

                for i in xrange(len(seq_a)):
                    if not seq_a[i] == seq_q[i]:
                        seq_a[i] = seq_a[i].lower()
                res[id_q].seq_q = seq_q
                tracks = res[id_q].setdefault("tracks",{})
                tracks[trackName] = seq_a.tostring()
        nameFldLen = 20
        out = open(self.crArrBlastTxtFile,'w')
        for id_q,rec in sorted(res.items()):
            if not (opt.virMicTracksUnion and opt.virArrHitsOnly and "01.vir" not in rec.tracks):
                seqNames = ["Array".ljust(nameFldLen) ]
                seqs = [ rec.seq_q ]
                for trackName in sorted(hspsTracks.keys()):
                    seqNames.append(trackName.ljust(nameFldLen))
                    if rec.tracks.has_key(trackName):
                        seqs.append(rec.tracks[trackName])
                    else:
                        seqs.append('')
                out.write("Array %s\n" % id_q)
                printAliSeqs(seqs=seqs,lineLen=100,out=out,seqNames=seqNames,emptySymb='-')
                out.write("\n\n")
                #out.write("%s\n" % (seq_a.tostring(),))
        out.close()

    def exportBlastPerMatch(self,seqsQ,hspsTracks):
        res = defdict(Struct)
        for (trackName,hsps) in sorted(hspsTracks.items()):
            hsps = groupRecArray(hsps,"id_q")
            for (id_q,hsps_q) in sorted(hsps.items()):
                seq_q = seqsQ[id_q]
                seq_a_t = n.zeros(len(seq_q),dtype='S1')
                seq_a_t[:] = '.'
                #out.write("\n%s\n" % seq_q)
                hsps_q = groupRecArray(hsps_q,"id_h")
                for (id_h,hsps_h) in sorted(hsps_q.items()):
                    seq_a = seq_a_t.copy()
                    for hsp in hsps_h:
                        seq_ali_h = n.fromstring(hsp["seq_ali_h"],dtype='S1')
                        try:
                            seq_a[hsp["begin_q"]:hsp["end_q"]] = seq_ali_h
                        except:
                            ali_max_len = fieldDtypeRecArray(hsps_h,"seq_ali_h").itemsize
                            ali_q_len = -hsp["begin_q"] + hsp["end_q"]
                            if ali_q_len > ali_max_len:
                                print "Alignment truncated due to limited field length"
                                seq_a[hsp["begin_q"]:hsp["begin_q"]+ali_max_len] = seq_ali_h
                            else:
                                print "Unexplained length mismatch in alignment, truncating: %s %s" % (ali_q_len,len(seq_ali_h))
                            seq_a[hsp["begin_q"]:hsp["begin_q"]+len(seq_ali_h)] = seq_ali_h

                    for i in xrange(len(seq_a)):
                        if not seq_a[i] == seq_q[i]:
                            seq_a[i] = seq_a[i].lower()
                    key = (id_q,id_h)
                    res[key].seq_q = seq_q
                    tracks = res[key].setdefault("tracks",{})
                    tracks[trackName] = seq_a.tostring()
        nameFldLen = 20
        out = open(self.crArrBlastPerMatchTxtFile,'w')
        for (id_q,id_h),rec in sorted(res.items()):
            seqNames = ["Array".ljust(nameFldLen) ]
            seqs = [ rec.seq_q ]
            for trackName in sorted(hspsTracks.keys()):
                seqNames.append(trackName.ljust(nameFldLen))
                if rec.tracks.has_key(trackName):
                    seqs.append(rec.tracks[trackName])
                else:
                    seqs.append('')
            out.write("Array %s\tMatch %s\n" % (id_q,id_h))
            printAliSeqs(seqs=seqs,lineLen=100,out=out,seqNames=seqNames,emptySymb='-')
            out.write("\n\n")
            #out.write("%s\n" % (seq_a.tostring(),))
        out.close()

    def exportBlastMatches(self,hsps):
        raise NotImplementedError()
        def _defline_fields(fields,rec):
            """Create a string for FASTA defline with array fields"""
            return ' '.join([ "/%s=%s".replace(' ','_') % (field,rec[field]) for field in fields ])
        hsps = groupRecArray(hsps,"id_h")
        out = open(self.crArrBlastMatchFile,'w')
        blastDb = BlastDb(blastDataDir=self.opt.readsBlastDir)
        reader = blastDb.fastaReader(dbName=self.readsBlastDb,giFile=None,defLineTargetOnly=True,maxDegen=None)
        for rec in reader.records():
            id_db = rec.header().split()[1]
            if id_db in hsps:
                print id_db
                hsps_h = hsps[id_db]
                descr = ",".join([ _defline_fields(("id_q","begin_q","end_q","begin_h","end_h","forward","mism"),hsp) \
                        for hsp in hsps_h ])
                id = id_db
                seq = Seq(rec.sequence(),alphabet="dna")
                seq = SeqRecord(seq,id=id,description=descr)
                SeqIO.write([seq],out,"fasta")
        out.close()
            #id_arr = int(rec.header()[1:])
            #seq = rec.sequence()
        #for (id_h,hsps_h) in sorted(hsps.items()):

    def loadArrayBlastBin(self,inFile,minAlignLen=20,minBitScore=40,maxEvalue=1000,maxMism=1,out=None,debug=False):
        """Load into memory numpy record array previously created by parseArrayBlast().
        The records will be filtered according to cutoff parameters."""
        recs = loadObj(inFile)
        return recs[logicalAnd(\
                recs["align_len"]>=minAlignLen,
                recs["bits"]>=minBitScore,
                recs["expect"]<=maxEvalue,
                recs["mism"]<=maxMism\
                )]


    def parseArrayBlast(self,inFile,minAlignLen=10,minBitScore=10,maxEvalue=100,maxMism=100,out=None,debug=False):
        """Parse BLAST XML output and return as numpy record array.
        The records will be filtered according to cutoff parameters."""
        dtype=[("id_q","i4"),
            ("id_h","S%s"%self.maxSeqIdLen),
            ("begin_q","i8"),
            ("end_q","i8"),
            ("begin_h","i8"),
            ("end_h","i8"),
            #("seq_ali_h","O"),
            #("seq_ali_m","O"),
            ("seq_ali_h","S%s"%self.maxArrElemLen),
            ("seq_ali_m","S%s"%self.maxArrElemLen),
            ("forward","b"),
            ("mism","i8"),
            ("align_len","i8"),
            ("bits","i8"),
            ("expect","f4")]
        if out is None:
            out = sys.stdout
        inp = openCompressed(inFile,"r")
        blastRecs = NCBIXML.parse(inp)
        spcnt = defdict(int)
        arr = []
        revCompl = RevCompl()
        iRec = 0
        for blastRec in blastRecs:
            for alignment in blastRec.alignments:
                for hsp in alignment.hsps:
                    mism = hsp.align_length - hsp.identities
                    if hsp.align_length >= minAlignLen \
                            and hsp.bits >= minBitScore \
                            and hsp.expect <= maxEvalue \
                            and mism <= maxMism:
                    #if hsp.align_length >= 20 and \
                    #        abs(hsp.align_length - blastRec.query_letters) <= 2 and \
                    #        hsp.align_length - hsp.identities <= 2 and \
                    #        hsp.bits >= 40:
                        #print str(sorted(hsp.__dict__.items()))
                        #print str(sorted(blastRec.__dict__.items()))
                        #print str(sorted(alignment.__dict__.items()))
                        spcnt[blastRec.query] += 1
                        if debug:
                            print >> out, 'Query: ', blastRec.query
                            print >> out, 'Hit: ', alignment.hit_def
                            print >> out, hsp.query
                            print >> out, hsp.match
                            print >> out, hsp.sbjct
                            print >> out,   'Query lenght: ', blastRec.query_letters,\
                                    'Alignment length: ', hsp.align_length,\
                                    'Bit-score: ', hsp.bits,\
                                    'Query Coverage: ', round(float(hsp.align_length)/blastRec.query_letters*100),\
                                    'Mismatches: ', (hsp.align_length - hsp.identities),\
                                    'Identity: ', round(float(hsp.identities)/hsp.align_length*100),\
                                    'E-value: %e' % hsp.expect
                            print >> out, "\n\n"
                        id_q = blastRec.query
                        id_h = alignment.hit_def.split()[0]
                        assert hsp.sbjct_start < hsp.sbjct_end,"Only expecting reverse compliment query hsp coords in BLAST output, not hit hsp coords"
                        assert hsp.query_start != hsp.query_end,"I cannot determine the rev-compl state of one-letter hsp yet"
                        if hsp.query_start < hsp.query_end:
                            begin_q = hsp.query_start - 1
                            end_q =  hsp.query_end
                            begin_h = hsp.sbjct_start - 1
                            end_h = hsp.sbjct_end
                            #BIO parser returns unicode strings for some reason,
                            #and then translate in revCompl does not work
                            seq_ali_h = str(hsp.sbjct)
                            seq_ali_m = str(hsp.match)
                            forward = True
                        else:
                            begin_q = hsp.query_end - 1
                            end_q =  hsp.query_start
                            begin_h = alignment.length - hsp.sbjct_end
                            end_h = alignment.length - hsp.sbjct_start + 1
                            seq_ali_h = revCompl(str(hsp.sbjct))
                            seq_ali_m = str(hsp.match[::-1])
                            forward = False
                        arr.append((id_q,id_h,begin_q,end_q,begin_h,end_h,
                            seq_ali_h,seq_ali_m,forward,mism,
                            hsp.align_length,hsp.bits,hsp.expect))
                        iRec += 1
        arr = n.asarray(arr,dtype=dtype)
        #arr = n.rec.fromrecords(arr,dtype=dtype)
        if debug:
            print >> out, "Number of unique query sequencies with conforming hits: %s" % len(spcnt)
            print >> out, sorted(spcnt.items())
        return arr


    def loadArraySeqFile(self):
        reader = FastaReader(self.crArrSeqFile)
        arr = {}
        for rec in reader.records():
            id_arr = int(rec.header()[1:])
            seq = rec.sequence()
            arr[id_arr] = seq
        return arr

    def listPilerIoFiles(self,input=True):
        opt = self.opt
        if input:
            pilerSfx = ".fasta-*"
            return ( f for f in iglob(pjoin(opt.crArrSeqDir,"*"+pilerSfx)) if re.match(r'.*-[0-9]+$',f) )
        else:
            pilerSfx = ".piler.csv"
            return iglob(pjoin(opt.crArrDir,"*"+pilerSfx))

    def importPiler(self):

        opt = self.opt
        
        FLD_IND_ID = 0
        FLD_IND_POS = 1
        FLD_IND_SPLEN = 4
        FLD_IND_REP = -2
        FLD_IND_SP = -1

        assert opt.minSpacerLen > 0

        dtype=[("seq_part","S%s"%self.maxSeqPartIdLen),
            ("id_seq","S%s"%self.maxSeqIdLen),
            ("id_arr","i4"),
            ("pos","i8"),
            ("repeat","O"),
            ("spacer","O")]
        data = []
        id_arr = 0
        id_seqLast = ""
        
        inFiles = self.listPilerIoFiles(input=False)
        pilerSfx = ".piler.csv"
        
        for inFile in inFiles:
            # unique input seq file partition  is everything after the last '-' and before the piler.csv suffix, as in
            # asm.scf-00001.piler.csv
            inpPartId = stripSfx(os.path.basename(inFile),pilerSfx).rsplit('-',1)[1].strip() 
            inp = openCompressed(inFile,'r')

            for line in inp:
                parts = [ part.strip() for part in line.split('\t') ]
                spLen = int(parts[FLD_IND_SPLEN])
                spLen = int(parts[FLD_IND_SPLEN])
                id_seq = parts[FLD_IND_ID].split()[0]
                rep = parts[FLD_IND_REP]
                sp = parts[FLD_IND_SP]
                if id_seq != id_seqLast:
                    id_seqLast = id_seq
                data.append((inpPartId,id_seq,id_arr,int(parts[FLD_IND_POS])-1,rep,sp))
                if spLen < 0:
                    id_arr += 1
            inp.close()

        data = n.asarray(data,dtype=dtype)
        arr = {}
        for iRec in xrange(len(data)):
            rec = data[iRec]
            key = (rec['id_seq'].item(),rec['id_arr'].item())
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
            if (len(val) - 1) < opt.minSpacerNum:
                print "Number of spacers is less than %s in array, deleting:\n%s\n%s" % (opt.minSpacerNum,key,val)
                del arr[key]
            elif crisprHasSimilarSpacers(val):
                print "Found putative array with similar spacers, deleting:\n%s\n%s" % (key,val)
                del arr[key]
            elif not crisprLengthRatioOk(val):
                print "Found putative array that violates length ratios, deleting:\n%s\n%s" % (key,val)
                del arr[key]
        return arr
    
    def pilerToSql(self,arr):
        db = self.getDbSql()
        self.createTablesCrispr()
        idGenArr = IntIdGenerator()
        idGenArrElem = IntIdGenerator()
        inserterArr = db.makeBulkInserterFile(table='arr',bufLen=50000,workDir=self.tmpDir)
        inserterArrElem = db.makeBulkInserterFile(table='arr_elem',bufLen=50000,workDir=self.tmpDir)
        self.loadPilerArrSql(db=db,idGenArr=idGenArr,inserterArr=inserterArr,
                idGenArrElem=idGenArrElem,inserterArrElem=inserterArrElem,arr=arr)
        inserterArr.flush()
        inserterArrElem.flush()
        db.createIndices(table="arr",
            names=["id_seq"],
            primary="id")
        db.ddl("analyze table arr",ifDialect="mysql")
        db.createIndices(table="arr_elem",
            names=["id_arr"],
            primary="id")
        db.ddl("analyze table arr_elem",ifDialect="mysql")
        self.delDbSql()

    def loadPilerArrSql(self,db,idGenArr,inserterArr,
            idGenArrElem,inserterArrElem,arr):
        """Load array dict prepared by importPiler() into SQL tables.
        """
        for ((id_seq,id_arr_old),recs) in sorted(arr.items()):
            id_arr = idGenArr()
            seq_part = recs[0]["seq_part"]
            inserterArr((id_arr,id_seq,seq_part))
            for (iRec,rec) in izipCount(recs):
                begin = rec['pos']
                seq_ali = rec['repeat']
                end = begin + len(crisprSeqAliToRaw(seq_ali))
                id_elem = idGenArrElem()
                inserterArrElem((id_elem,id_arr,begin,end,1,seq_ali))
                if iRec < len(recs) - 1:
                    begin = end
                    seq_ali = rec['spacer']
                    end = begin + len(crisprSeqAliToRaw(seq_ali))
                    id_elem = idGenArrElem()
                    inserterArrElem((id_elem,id_arr,begin,end,0,seq_ali))

    def loadCrisprSqlArr(self):
        """Load CRISPR from SQL tables into an array"""
        db = self.getDbSql()
        recs = db.selectAsArray("""
        select *
        from arr_elem
        order by id_arr,begin
        """)
        self.delDbSql()
        return recs


    def createTablesCrispr(self):
        db = self.getDbSql()
        db.ddl("""
        create table arr
        (
        id integer,
        id_seq char(%s),
        seq_part char(%s)
        )
        """ % (self.maxSeqIdLen,self.maxSeqPartIdLen),
        dropList=["table arr"])
        db.ddl("""
        create table arr_elem
        (
        id integer,
        id_arr integer,
        begin integer,
        end integer,
        is_rep bool,
        seq_ali varchar(%s)
        )
        """ % (self.maxArrElemLen,),
        dropList=["table arr_elem"])

    def getTaxaTree(self):
        if self.taxaTree is None:
            self.taxaTree = loadTaxaTree()
        return self.taxaTree

    def exportCrisprArraySeq(self):
        db = self.getDbSql()
        recs = db.selectAsArray("""
        select *
        from arr_elem
        order by id_arr,begin
        """)
        self.delDbSql()
        recs = groupRecArray(recs,"id_arr")
        out = open(self.crArrSeqFile,'w')
        for (id_arr,recs_arr) in sorted(recs.items()):
            seq = []
            for rec in recs_arr:
                s = crisprSeqAliToRaw(rec["seq_ali"])
                if rec["is_rep"]:
                    s = s.lower()
                else:
                    s = s.upper()
                seq.append(s)
            seq = Seq(''.join(seq),alphabet="dna")
            seq = SeqRecord(seq,id="%s" % (id_arr,),description='')
            SeqIO.write([seq],out,"fasta")
        out.close()



    def loadCrArrayRangesSql(self,db):
        """Return numpy record array from SQL tables with total coord range of every CRISPR array"""
        return db.selectAsArray("""
        select a.id_seq,b.id_arr,min(b.begin) as begin,max(b.end) as end 
        from arr a, arr_elem b 
        where a.id=b.id_arr group by a.id_seq,b.id_arr
        """)

    def loadCrArrayElemSql(self,db):
        """Return numpy record array from SQL tables with every element of every CRISPR array.
        id_seq is included, and records are sorted by (id_arr,begin)."""
        return db.selectAsArray("""
        select a.id_seq,b.id_arr,b.id as id_elem,b.begin,b.end,b.is_rep
        from arr a, arr_elem b 
        where a.id=b.id_arr
        order by b.id_arr,b.begin
        """)
    
    def getStatFileName(self,sfx):
        return self.statFileRoot+'.'+sfx+'.csv'

    def stats(self):
        db = self.getDbSql()
        db.exportAsCsv("select * from st_org_all","tmp_st_org_all.csv")
        self.delDbSql()
        return
        # All arrays with total span
        db.createTableAs("st_arr",
        """
        select a.id_seq,b.id_arr,min(b.begin) as begin, max(b.end) as end,(count(*)+1)/2 as cnt_rep
        from arr a, arr_elem b
        where a.id = b.id_arr
        group by a.id_seq,b.id_arr
        """)
        # All sequences with array counts
        db.createTableAs("st_seq_arr",
        """
        select id_seq,count(*) as cnt_arr,max(cnt_rep) as max_cnt_rep
        from st_arr
        group by id_seq
        """)
        # All sequences with min distance between CAS and CRISPR
        # A distance between midpoints of a protein and an array minus a sum of their half-length
        # is the min distance between their ends. We compute this distance so that later we could
        # impose a cutoff of, say, 500 bp and see which sequences have CAS proteins and CRISPR array 
        # located close to each other.
        db.createTableAs("st_prot_arr",
        """
        select a.id_seq,min(abs((a.begin  + a.end)/2 - (b.begin+b.end)/2) - abs(a.end-a.begin)/2 - abs(b.end-b.begin)/2) as prot_arr_dist
        from mic_prot_annot a, st_arr b
        where a.id_seq = b.id_seq
        group by a.id_seq
        """)
        # All sequences with CAS counts
        db.createTableAs("st_seq_prot",
        """
        select id_seq,count(*) as cnt_prot
        from mic_prot_annot
        group by id_seq
        """)
        # All sequences left joined with CRISPR counts, CAS counts and min CRISPR-CAS distances
        db.createTableAs("st_seq_all",
        """
        select a.name,a.id_seq,b.cnt_arr,b.max_cnt_rep,c.cnt_prot,d.prot_arr_dist,a.src,a.genel,a.seq_len
        from mic_seq a
        left join st_seq_arr b
        on a.id_seq = b.id_seq
        left join st_seq_prot c
        on a.id_seq = c.id_seq
        left join st_prot_arr d
        on a.id_seq = d.id_seq
        """)
        # All organisms left joined with CRISPR counts, CAS counts and min CRISPR-CAS distances
        db.createTableAs("st_org_all",
        """
        select name,src,sum(cnt_arr) as cnt_arr,
        max(max_cnt_rep) as max_cnt_rep,
        sum(cnt_prot) as cnt_prot,
        min(prot_arr_dist) as prot_arr_dist,
        sum(seq_len) as sum_seq_len,
        round(avg(seq_len)) as avg_seq_len,
        count(*) as cnt_seq
        from st_seq_all 
        group by name,src
        """)


def run_Crispr():
    opt = Struct()
    opt.runMode = "inproc" #"batchDep"
    #modes = [ "blastcr","loadbl","annotmic","stats" ] 
    modes = [ "annotmic" ]
    #modes = [ "loadmic","blastdb","pilercr","loadcr","exportcr","blastcr","loadbl","annotmic","stats" ] 
    #modes = [ "stats" ] #"annotmic" "loadmicprot" "annotmic" "pilercr" "loadcr" "exportcr" "loadmic" "loadannot" "exportaclame" "exportmgt" "blastcr" "loadbl" "blastdb"
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = CrisprApp(opt=opt)
        print "Running mode %s" % mode
        jobs = app.run(depend=jobs)

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(CrisprApp)
    #run_Crispr()

