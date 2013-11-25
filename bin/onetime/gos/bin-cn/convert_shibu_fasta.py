from MGT.Util import *
from MGT.FastaIO import *
from MGT.Taxa import *
from MGT.UUID import *
import json, glob


models = []

fastaWriter = FastaWriter(sys.argv[2])

taxaTree = loadTaxaTree()

for inp in glob.glob(sys.argv[1]):

    print "Reading {}".format(inp)

    fastaReader = FastaReader(inp)

    for rec in fastaReader.records():
        hdr = rec.header()
        seqid = rec.getId() # TAXIDXXXX_YYYYY
        assert seqid.lower().startswith("taxid")
        taxid = int(seqid[len("taxid"):].split("_",1)[0])
        try:
            node = taxaTree.getNode(taxid)
        except:
            print "taxid {0} not found".format(seqid)
            raise
        name = seqid.split("_",1)[1]
        seqid = name[:maxIdLen]
        hdr = seqid
        fastaWriter.record(header=hdr,sequence=rec.sequence())
        model = dict(
                taxid = taxid,
                id = seqid,
                name = name,
                ids_seq = [ seqid ]
                )
        models.append(model)

fastaWriter.close()
with open(sys.argv[3],"w") as modelsOut:
    json.dump(models,modelsOut)


