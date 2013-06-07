from MgxAnnotFormatApisCommon import *
from MgxAnnotFormatCommon import *

_apisAnnotRecSamp = """\
Seq Name        Dataset Length  Annotation      Standard Classification Standard Tax Id
scf7180008708094_1393_1809_2    CN_GOS4_92      139     homocysteine_S-methyltransferase_family_protein Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Marinobacter; Mixed   2742
scf7180008708094_445_1_1        CN_GOS4_92      148             NO_TREE
scf7180008708095_10079_9078_12  CN_GOS4_92      334     hypothetical_protein    Bacteria; Proteobacteria; Mixed 1224
scf7180008708095_10393_10076_13 CN_GOS4_92      106     transcriptional_regulator_BolA  Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Saccharophagus; Saccharophagus degradans      86304
scf7180008708095_1060_3723_2    CN_GOS4_92      888     glnD    Bacteria; Proteobacteria; Gammaproteobacteria; Mixed    1236
scf7180008708095_10892_10415_14 CN_GOS4_92      159     dTDP-glucose_4,6-dehydratase_(EC:4.2.1.46)      Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas; Pseudomonas fluorescens  294
scf7180008708095_275_1060_1     CN_GOS4_92      262     map     Mixed
scf7180008721125_1252_308_2     CN_GOS4_92      315     iunH    Eukaryota; Bacillariophyta; Bacillariophyceae; Bacillariaceae (order); Bacillariaceae; Nitzschia; Nitzschia I-146       -283592
"""

class ApisAnnotRecord(object):

    def __init__(self,rec):
        self.rec = rec
        self.parsed = {"id_q":rec[0].strip()}

    def parse(self):
        parsed = self.parsed
        rec = self.rec
        try:
            parsed["taxid"] = int(rec[-1].strip())
        except:
            pass
        else:
            parsed["id_s"] = ""
            parsed["descr"] = rec[3].strip()
            parsed["conf_x"] = 1.
            tax_lin_names = [ w.strip() for w in rec[4].strip().split(";") ]
            parsed["tax_lin_names"] = tax_lin_names
            parsed["tax_name"] = tax_lin_names[-1]
            for tax_name in reversed(tax_lin_names):
                if tax_name.lower() != "mixed":
                    parsed["tax_name"] = tax_name
                    break
        return parsed

class ApisAnnotReader(object):
    
    def __init__(self,inp,pepMap=None):
        x = openCsv(inp,
                mode="r",
                factory=csv.reader,
                dialect="excel-tab")
        self.reader = x.csvFile
        self.csvClose = x.csvClose
        self.csvFileInp = x.csvFileStream
        self.pepMap = pepMap

    def next(self):
        row = self.reader.next()
        #skip header
        while row[0].strip() == "Seq Name":
            row = self.reader.next()
        rec = ApisAnnotRecord(row)
        rec.parse()
        if self.pepMap:
            self._updateFromPep(rec)
        else:
            rec.parsed.update(parseOrfOnContigIdApis(rec.parsed["id_q"]))
        return rec
    
    def __iter__(self):
        return self

    def close(self):
        if self.csvClose:
            self.csvFileInp.close()
           
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
