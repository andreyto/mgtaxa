### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Utilities to query NCBI Entrez through Eutils.
Built on top of Biopython."""

from MGT.Common import *

from Bio import Entrez, SeqIO

def _idJoin(ids):
    return ",".join(["%s" % id for id in ids])

class EzRequest:
    """Perform large batch requests to Entrez.
    Post a list of IDs, 
    download records in batches,
    cache results in a file,
    parse and iterate through results.
    """

    def __init__(self,toolName=options.toolName,
            email=options.toolEmail,
            batchSize=1000,
            cacheDir=None,
            db="nucleotide",
            rettype="genbank"):
        #rettype="genbak" produces deprecation warning from Biopython,
        #"gb" or "gp" that it advises return nothing.
        self.toolName = toolName
        self.email = email
        self.batchSize = batchSize
        if not cacheDir:
            cacheDir = os.getcwd()
        self.cacheDir = cacheDir
        self.db = db
        self.rettype = rettype
        makedir(cacheDir)

    def post(self,ids):
        print "DEBUG: Posting the GI list to Entrez session"
        res = Entrez.read(\
                Entrez.epost(db=self.db, id=_idJoin(ids),tool=self.toolName,email=self.email)\
                )
        return res


    def fetch_blocked(self,ids,outFile=None,seq_range=(1,1),**kw):
        """Fetch records from Entrez.
        @param ids list of primary IDs (GIs for sequence records). 
        @param outFile Output file name or file-like object
        @param seq_range How much to fetch from each sequence.

        Eutils have no post function for batch use of ACC.
        Also, epost used here does not recognize aggregate "sequences" name for DB.
        """
        res = self.post(ids)
        webenv = res["WebEnv"]
        query_key = res["QueryKey"]
        ownFile = False
        if outFile is None:
            out_handle,outFile = makeTmpFile(".req",dir=self.cacheDir)
            ownOut = True
            ownFile = True
        elif not hasattr(outFile,"write"):
            out_handle = open(outFile, "w")
            ownOut = True
        else:
            out_handle = outFile
            ownOut = False
        
        for retstart in range(0,len(ids),self.batchSize):
            retstart = 0
            net_handle = Entrez.efetch(db=self.db,
                    rettype=self.rettype,
                    retstart=retstart,
                    retmax=self.batchSize,
                    webenv=webenv, 
                    query_key=query_key,
                    tool=self.toolName,
                    email=self.email,
                    seq_start=seq_range[0],
                    seq_stop=seq_range[1],
                    **kw)
            buff_size = 2**16
            while True:
                buff = net_handle.read(buff_size)
                if not buff:
                    break
                out_handle.write(buff)
            net_handle.close()
        if ownOut:
            out_handle.close()
        return outFile


    def fetch_blocked_parsed(self,ids,outFile=None,seq_range=(1,1),**kw):
        assert outFile is None or not hasattr(outFile,"write"),\
                "Cannot parse when saving to existing output stream"
        outF = self.fetch_blocked(ids=ids,outFile=outFile,seq_range=seq_range,**kw)
        in_handle = open(outF,'r')
        if self.rettype == "xml":
            rec = Entrez.read(in_handle)
            yield rec
            #XML return set is not a list, parse() raises an exception.
            #for rec in Entrez.read(in_handle):
            #    yield rec
        elif self.rettype == "genbank":
            for rec in SeqIO.parse(in_handle,"genbank"):
                yield rec
        in_handle.close()
        if outFile is None:
            os.remove(outF)

    def summary(self,ids,**kw):
        res = self.post(ids)
        webenv = res["WebEnv"]
        query_key = res["QueryKey"]
        net_handle = Entrez.esummary(db=self.db,
                webenv=webenv, 
                query_key=query_key,
                tool=self.toolName,
                email=self.email,
                **kw)
        for rec in Entrez.parse(net_handle):
            yield rec
        net_handle.close()
