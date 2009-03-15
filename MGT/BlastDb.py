"""Interface into NCBI BLAST DB: allows to extract FASTA sequence data for a given database or a list of GIs.
Calls NCBI 'fastacmd' executable.
"""

from MGT.Common import *
from MGT.FastaIO import *

class BlastDb(MGTOptions):

    def __init__(self):
        MGTOptions.__init__(self)
        self.ncbiDbs =  (
                         Struct(id='g',db='refseq_genomic'),
                         Struct(id='o',db='other_genomic'),
                         Struct(id='w',db='wgs'),
                         Struct(id='h',db='htgs'),
                         Struct(id='n',db='nt')
                         )
        self.taxaTree = None

    def getDbs(self):
        return self.ncbiDbs

    def giFileToBin(self,giFile,giFileBin):
        giFile = os.path.abspath(giFile)
        run(("formatdb -F %s -B %s" % (giFile,giFileBin)).split(),cwd=self.blastDataDir)

    def makeDbAlias(self,dbNameAlias,giFile=None,dbNames=None):
        from glob import glob
        from datetime import datetime
        from textwrap import dedent
        if dbNames is None:
            dbNames = [ db.db for db in self.ncbiDbs ]
        directDbRefs = False
        giFileBin = dbNameAlias+'.gil'
        #comment out GILIST field in alias file if giFile is None
        if giFile is None:
            comGILIST = '#'
        else:
            comGILIST = ''
            self.giFileToBin(giFile,giFileBin)
        aliasStr = dedent("""\
        #
        # Alias file created %s
        #
        #
        TITLE MGT sequence selection
        #
        DBLIST %%s
        #
        %sGILIST %s
        #
        #OgiFile
        #
        """ % (datetime.today().ctime(),comGILIST,giFileBin))
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

    def fastaStream(self,dbName,giFile=None,defLineTargetOnly=True,maxDegen=20):
        """Return an open file object that streams database records in FASTA format.
        @param dbName - database name (e.g. 'nt')
        @param giFile - text file with gi ids (one per line)
        @defLineTargetOnly - switch on '-t' option of fastacmd. Some records have
        several different deflines pointing to the same sequence (most such records are the result
        of creating RefSeq entries from unmodified GenBank entries, but not all). Without this switch set,
        'fastacmd' will output a defline as in the following example, for either 23455713 or 15183 as a query gi:
        >gi|23455713|ref|NC_004301.1| Enterobacteria phage FI, complete genome >gi|15183|emb|X07489.1| Bacteriophage SP genomic RNA
        For each query gi, only one record will be output. If this switch is set, then defline will only contain information
        related for the query gi, as in:
        >gi|23455713|ref|NC_004301.1| Enterobacteria phage FI, complete genome
        Note: 'fastacmd' will refuse to dump (-D 1) database whose alias file includes gi list. The only option
        is to use '-i' instead and provide a text gi list file"""
        otherOpts = ""
        if defLineTargetOnly:
            otherOpts = otherOpts + "-t T"
        if giFile is None:
            giListOpt = "-D 1"
        else:
            giFile = os.path.abspath(giFile)
            giListOpt = "-i %s" % (giFile,)
        # -l 60000 seqfaults fastacmd 2.2.15 SuSe 10 x86_64 (boxer)
        if self.debugFakeSequence:
            fastaCmdBin = "fastacmd_debug"
        else:
            fastaCmdBin = "fastacmd"
        cmd = fastaCmdBin + " -l 30000 -d %s %s %s" % (dbName,giListOpt,otherOpts)
        # Initially I implemented the pipe using 'Popen shell pipe' idiom (see subprocess module help),
        # but it was x10 slower than the true shell pipe below.
        if maxDegen is not None:
                cmd = cmd + " | " + os.path.join(os.environ["MGT_EXEC_BIN"],"fastaFilter") + " -n %s -l %s" % (maxDegen,2**16)
        if self.debug > 0:
            print cmd
        return Popen( cmd,
                        cwd=self.blastDataDir,
                        env=os.environ,
                        bufsize=2**16,
                        stdout=PIPE,
                        shell=True,
                        close_fds=True).stdout

    def fastaReader(self,**kw):
        return FastaReader(self.fastaStream(**kw))

