from MgxAnnotFormatApisCommon import *
from MgxAnnotFormatCommon import *
import MGT.JSON

_apisAnnotRecSamp = """\
{
  "name": "JCVI_PEP_metagenomic.orf.1134269553525.1",
  "tree": "((AeJTIv/++Di/jskyI9rCelRGskY__UDEL_CHESAPEAKE_VIROPLANKTON_SMPL:0.30162,
  (L32WHWBUuzwH975IX6EAYl1hP48__UDEL_CHESAPEAKE_VIROPLANKTON_SMPL:0.25973,
  h9X4Sj60iJzUlCSchszbOLozkqk__Psychroflexus_torquis_ATCC_700755:0.34433)69:0.0235)89:0.03972,
  (JCVI_PEP_metagenomic.orf.1134269553525.1:0.28773,
  (((GqPn4m8LUeR97mPo2qCZNAkN3LE__Burkholderia_sp._Ch1-1,GqPn4m8LUeR97mPo2qCZNAkN3LE__Burkholderia_sp._Ch1-1):
  0.24377,scQyBL12xTPWM2XZ80FHueZcBUQ__Desulfovibrio_sp._FW1012B:0.29177)74:0.03626,
  (PgflL12p8FsipcqlgjRC3rRXBvo__UDEL_CHESAPEAKE_VIROPLANKTON_SMPL:0.04049,
  (qHWBb9wzVp2Qm6Lhhe5/627zIVg__UDEL_CHESAPEAKE_VIROPLANKTON_SMPL,
  qHWBb9wzVp2Qm6Lhhe5/627zIVg__UDEL_CHESAPEAKE_VIROPLANKTON_SMPL):0.05326)100:0.24203)
  :0.00388)33:0.00379999999999997);",
  "strict_consensus": {
    "kingdom": "Bacteria",
    "kingdom_outgroup": 0,
    "phylum": "Proteobacteria",
    "phylum_outgroup": 1,
    "class": "Mixed",
    "class_outgroup": 0,
    "order": "Mixed",
    "order_outgroup": 0,
    "family": "Mixed",
    "family_outgroup": 0,
    "genus": "Mixed",
    "genus_outgroup": 0,
    "species": "Mixed",
    "species_outgroup": 0
  },
  "relaxed_consensus": {
    "kingdom": "Bacteria",
    "kingdom_outgroup": 0,
    "phylum": "Proteobacteria",
    "phylum_outgroup": 1,
    "class": "Betaproteobacteria",
    "class_outgroup": 1,
    "order": "Burkholderiales",
    "order_outgroup": 1,
    "family": "Burkholderiaceae",
    "family_outgroup": 1,
    "genus": "Burkholderia",
    "genus_outgroup": 1,
    "species": "Burkholderia sp. Ch1-1",
    "species_outgroup": 1
  },
  "strict_taxon_id": "1224",
  "strict_taxon": "Proteobacteria",
  "relaxed_taxon_id": "243261",
  "relaxed_taxon": "Burkholderia sp. Ch1-1",
  "bootstrap": "33"
}
"""

class ApisAnnotRecord(object):
    
    lin_ranks = ("kingdom","phylum","class","order","family","genus","species")

    taxon_call_prefix = "relaxed"

    def __init__(self,rec):
        self.rec = rec
        self.parsed = {"id_q":rec["name"]}

    def parse(self):
        parsed = self.parsed
        rec = self.rec
        try:
            parsed["taxid"] = int(rec[self.taxon_call_prefix + "_taxon_id"])
        except:
            pass
        else:
            parsed["id_s"] = ""
            parsed["descr"] = rec.get("annotation","")
            parsed["conf_x"] = 1.
            lin = rec[self.taxon_call_prefix + "_consensus"]
            lin_ranks = self.lin_ranks
            tax_lin_names = [lin[rank] for rank in lin_ranks if rank in lin and lin[rank].lower() != "mixed"]
            parsed["tax_lin_names"] = tax_lin_names
            parsed["tax_name"] = rec[self.taxon_call_prefix + "_taxon"]
        return parsed


class ApisAnnotReader(object):
    
    @staticmethod
    def sortInPepOrder(apisAnnotInp,apisAnnotOut,inpPep=None,dropTreeAttrib=True):
        import json
        assert inpPep is not None, "Currently inpPep must be always provided"
        inpPep = MappedPepFastaReader(inpPep)
        idPepOrd = dict( ( (rec.getId(),i) for (i,rec) in enumerate(inpPep.records())) )
        inpPep.close()
        #we really only need the reader attribute from apisAnnotInp
        apisInpStream = openCompressed(apisAnnotInp,"r")
        reader = MGT.JSON.json_records_iterator(apisInpStream,withJson=True)
        ordApis = []
        for i,(row,rowStr) in enumerate(reader): 
            if not row["tree"] == "NO_TREE":
                if dropTreeAttrib:
                    del row["tree"]
                    jsonStr = json.dumps(row)
                else:
                    jsonStr = rowStr
                ordApis.append((idPepOrd[row["name"]],jsonStr))
            if i % 10**5 == 0:
                print "DEBUG: loaded JSON APIS record with original index {0}".format(i)
        apisInpStream.close()
        ordApis.sort()
        writer = openCompressed(apisAnnotOut,"w")
        for (i,row) in ordApis:
            writer.write("{0}\n".format(row))
            if i % 10**5 == 0:
                print "DEBUG: wrote JSON APIS record with sorted index {0}".format(i)
        writer.close()
    
    def __init__(self,inp,pepMap=None):
        if isinstance(inp,str):
            inp = openCompressed(inp,"r")
            self.inpClose = True
        else:
            self.inpClose = False
        self.inp = inp
        self.reader = MGT.JSON.json_records_iterator(inp)
        self.pepMap = pepMap

    def next(self):
        for row in self.reader:
            #APIS JSON export puts "NO_TREE" records
            #all at the end, violating the order of the
            #peptide records. We should skip "NO_TREE"
            #records because they are empty anyway
            if not row.get("tree","") == "NO_TREE":
                #tree is a large field, wack it for now
                if "tree" in row:
                    del row["tree"]
                rec = ApisAnnotRecord(row)
                rec.parse()
                if self.pepMap:
                    self._updateFromPep(rec)
                else:
                    rec.parsed.update(parseOrfOnContigIdApis(rec.parsed["id_q"]))
                return rec
        else:
            raise StopIteration
    
    def __iter__(self):
        return self

    def close(self):
        if self.inpClose:
            self.inp.close()
           
    def _updateFromPep(self,rec):
        id_q = rec.parsed["id_q"]
        pepRec = self.pepMap.findIdPep(id_q)
        assert pepRec,"Peptide record not found in peptide map for peptide id: %s" % (id_q,)
        parsed = rec.parsed
        parsed["id_cont"] = pepRec["read_id"]
        parsed["start_orf"] = pepRec["begin"]
        parsed["end_orf"] = pepRec["end"]
        parsed["strand_orf"] = pepRec["orientation"]
        parsed["pep_id"] = pepRec["pep_id"]

ApisAnnotReaderPepIds = ApisAnnotReader

