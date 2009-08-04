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
from Bio import pairwise2
from Bio.Blast import NCBIXML

from glob import iglob


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
    """App-derived class for CRISPR array processing"""

    batchDepModes = ("blastcr",)

    maxSeqIdLen = UUID.maxIdLen
    maxSeqPartIdLen = UUID.maxIdLen
    maxArrElemLen = 200

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        opt.sqlHost = "atovtchi1-lx.jcvi.org"
        opt.minSpacerLen = 20
        opt.maxSpacerLen = 60
        opt.minSpacerNum = 3
        workDir = os.environ["GOS_WORK"]
        opt.topWorkDir = workDir
        opt.crDir = pjoin(workDir,"scaf-crispr")
        opt.crSeqDir = pjoin(workDir,"scaf-piler")
        opt.gosBlastDbDir = pjoin(workDir,"reads-blast")
        opt.crBlastDir = pjoin(workDir,"crispr-blast")
        opt.MEM = 6000
   
    def getDbSql(self):
        """Allocate (if necessary) and return a connection to SQL server"""
        if not hasattr(self,"dbSql"):
            self.dbSql = DbSqlMy(db="crispr",host=self.opt.sqlHost)
        return self.dbSql

    def delDbSql(self):
        """Call this to free a connection to SQL server if it will not be needed for extended period of time"""
        if hasattr(self,"dbSql"):
            self.dbSql.close()
            del self.dbSql

    def init(self):
        opt = self.opt
        self.taxaTree = None #will be lazy-loaded
        self.store = SampStore.open(path=self.opt.get("cwd",os.getcwd()))
        self.tmpDir = self.store.getFilePath("tmp")
        makedir(self.tmpDir)
        self.crArrSeqFile = self.store.getFilePath("arr.fasta")
        self.crArrBlastFileRoot = self.store.getFilePath("arr.blast")
        self.crArrBlastTxtFile = self.store.getFilePath("arr.blast.txt")
        self.crArrBlastPerMatchTxtFile = self.store.getFilePath("arr.blast.match.txt")
        self.crArrBlastMatchFile = self.store.getFilePath("arr.blast.match.fasta")
        self.initBlastJobs()
        self.mgtOutFile = "/usr/local/projects/GOSII/MGTAXA/v0.90/asm11_Mar06/genus/idlin.csv"
        self.aclameBlastFile = pjoin(opt.topWorkDir,"shannon/arr.blast.match_ACLAME.txt")
        self.aclameBlastOutFile = pjoin(opt.topWorkDir,"shannon/arr.blast.match_ACLAME.cov.txt")
   
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

    def doWork(self,**kw):
        self.init()
        opt = self.opt
        if opt.mode == "loadcr":
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
        elif opt.mode == "exportbl":
            return self.exportBlast()
        elif opt.mode == "exportmgt":
            return self.exportMgtaxa()
        elif opt.mode == "exportaclame":
            return self.exportAclameBlast()
        elif opt.mode == "loadannot":
            return self.loadAnnot()
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def loadCrispr(self):
        """Import and filter original PILERCR output files and convert them into SQL tables"""
        pass
        arr = self.importPiler()
        #self.store.saveObj(arr,"crispr_dict")
        #arr = self.store.loadObj("crispr_dict")
        self.pilerToSql(arr)

    def exportCr(self):
        self.exportCrisprArraySeq()

    def makeBlastDb(self):
        opt = self.opt
        curdir = os.getcwd()
        os.chdir(opt.gosBlastDbDir)
        # formatdb take a single space separatated string for -i option
        # when you need to combine several FASTA files to a single DB,
        # and it is crazily sensitive to delimiteres, that is why find -printf is used below
        # if we need to search in subdirectories.
        #seqFiles=$(find GSIO* -name inp.fasta -printf " %p")
        seqFiles=backsticks("""find . -name '*.fasta' -printf " %p" """,shell=True)
        gosBlastDb = os.path.basename([ db.dbPath for db in self.blastDbs if db.kind == "gos" ][0])
        run("""formatdb -n %s -t %s -p F -i "%s" """ % (gosBlastDb,gosBlastDb,seqFiles,),shell=True)
        os.chdir(curdir)

    def blastCr(self,**kw):
        opt = self.opt
        #pandaDbs = ["Phage.fasta"]
        jobs = []
        for blDb in self.blastDbs:
            blOpt = copy(opt)
            blOpt.mode = "blastcr-one"
            blOpt.blDbPath = blDb.dbPath
            blOpt.blOutFile = blDb.outFile
            blOpt.blInpFile = blDb.inpFile
            blOpt.MEM = 6000
            blOpt.LENGTH = "medium"
            blOpt.ARCH = "lx*64"
            blApp = self.factory(opt=blOpt)
            print blApp.opt
            jobs.append(blApp.run(**kw))
        return jobs

    def blastCrOne(self):
        opt = self.opt
        runBlast(dbPath=opt.blDbPath,inpFile=opt.blInpFile,outFile=opt.blOutFile,paramStr="-p blastn -e 1 -W 7 -F F -m 7 -U F") 

    def blastToCoverage(self):
        """Calculate how much spacers vs repeats each BLAST match hits.
        @return record array with (id_arr,id_hit,coverage_spacers,coverage_repeats)"""
        hsps = []
        for blDb in self.blastDbs:
            if blDb.group in ("vir",):
                hspsOne = self.parseArrayBlast(blDb.outFile,minAlignLen=22,minBitScore=1,maxEvalue=1000,maxMism=2,out=None)
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
        hspsTracks = {}
        for blDb in self.blastDbs:
            if blDb.group not in ("vir",):
                continue
            hsps = self.parseArrayBlast(blDb.outFile,minAlignLen=22,minBitScore=1,maxEvalue=1000,maxMism=0,out=None)
            #if blDb.group in ("ViralGroupNucl","Phage", "vir"):
            #    trackName  = "01.vir"
            #else:
            #    trackName = "10.mic"
            trackName = ("%2i"%blDb.order)+'.'+blDb.kind+'.'+blDb.group
            if hspsTracks.has_key(trackName):
                hspsTracks[trackName] = n.concatenate([hspsTracks[trackName],hsps])
            else:
                hspsTracks[trackName] = hsps
        seqsQ = self.loadArraySeqFile()
        #self.exportBlast(seqsQ=seqsQ,hspsTracks=hspsTracks)
        self.exportBlastPerMatch(seqsQ=seqsQ,hspsTracks=hspsTracks)
        #self.exportBlastMatches(hsps=hsps)

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

    def exportBlast(self,seqsQ,hspsTracks):
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
                            seq_a[hsp["begin_q"]:hsp["begin_q"]+len(seq_ali_h)] = seq_ali_h

                for i in xrange(len(seq_a)):
                    if not seq_a[i] == seq_q[i]:
                        seq_a[i] = seq_a[i].lower()
                res[id_q].seq_q = seq_q
                tracks = res[id_q].setdefault("tracks",{})
                tracks[trackName] = seq_a.tostring()
        nameFldLen = 20
        out = open(self.crArrBlastTxtFile,'w')
        for id_q,rec in sorted(res.items()):
            if not rec.tracks.has_key("01.vir"):
                continue
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

    def parseArrayBlast(self,inFile,minAlignLen=20,minBitScore=40,maxEvalue=1,maxMism=1,out=None,debug=False):
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
            ("mism","i4")]
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
                        arr.append((id_q,id_h,begin_q,end_q,begin_h,end_h,seq_ali_h,seq_ali_m,forward,mism))
                        iRec += 1
        arr = n.asarray(arr,dtype=dtype)
        #arr = n.rec.fromrecords(arr,dtype=dtype)
        if debug:
            print >> out, "Number of unique spacers with conforming hits: %s" % len(spcnt)
            print >> out, sorted(spcnt.items())
        return arr

    def loadAnnot(self):
        #self.loadPepToRead()
        #self.loadAsmToRead()
        #self.loadPepAnnot()
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

    def loadArraySeqFile(self):
        reader = FastaReader(self.crArrSeqFile)
        arr = {}
        for rec in reader.records():
            id_arr = int(rec.header()[1:])
            seq = rec.sequence()
            arr[id_arr] = seq
        return arr

    def importPiler(self):
        opt = self.opt
        pilerSfx = ".piler.csv"
        
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
        
        inFiles = iglob(pjoin(opt.crDir,"*"+pilerSfx))
        
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


def run_Crispr():
    opt = Struct()
    opt.runMode = "inproc" #"batchDep"
    modes = ["loadannot"] #"exportaclame" "exportmgt" "blastcr" "loadbl" "blastdb" "exportcr" "loadcr"
    jobs = []
    for mode in modes:
        opt.mode = mode
        app = CrisprApp(opt=opt)
        jobs = app.run(depend=jobs)

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(CrisprApp)
    #run_Crispr()

