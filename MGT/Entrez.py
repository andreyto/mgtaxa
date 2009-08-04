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
            cacheDir=os.path.join(options.tmpDir,"entrez"),
            db="nucleotide",
            rettype="genbank"):
        self.toolName = toolName
        self.email = email
        self.batchSize = batchSize
        self.cacheDir = cacheDir
        self.db = db
        self.rettype = rettype
        makedir(cacheDir)


    def fetch(self,ids,outFile=None,seq_range=(1,1)):
        """Fetch records from Entrez.
        @param list of primary IDs (GIs for sequence records). 
        Eutils have no post function for batch use of ACC.
        Also, epost used here does not recognize aggregate "sequences" name for DB.
        """
        print "Posting the GI list to Entrez session"
        res = Entrez.read(\
                Entrez.epost(db=self.db, id=_idJoin(ids),tool=self.toolName,email=self.email)\
                )
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
                    seq_stop=seq_range[1])
            out_handle.write(net_handle.read())
            net_handle.close()
            print "Fetched %s records out of %s" % (retstart+self.batchSize,len(ids))
        if ownOut:
            out_handle.close()
        return outFile


    def fetchParsed(self,ids,outFile=None,seq_range=(1,1)):
        assert outFile is None or not hasattr(outFile,"write"),\
                "Cannot parse when saving to existing output stream"
        outF = self.fetch(ids=ids,outFile=outFile,seq_range=seq_range)
        in_handle = open(outF,'r')
        if self.rettype == "xml":
            for rec in Entrez.read(in_handle):
                yield rec
        elif self.rettype == "genbank":
            for rec in SeqIO.parse(in_handle,"genbank"):
                yield rec
        in_handle.close()
        if outFile is None:
            os.remove(outF)

