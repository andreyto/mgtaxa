from MGT.Common import *
from MGT.Logging import *

import itertools, csv
from argh import ArghParser,arg

log = logging.getLogger(os.path.basename(sys.argv[0]))

def annot_line_to_rec(line):
    words = line.rstrip("\n").split("\t")
    return dict(itertools.izip(words[::2],words[1::2]))

def as_tab(annot_file,tab_file,recType="annotPred",lenSamp_min=0,cntHits_min=0):
    """Export specific record type as tabular file"""
    with closing(openCompressed(annot_file,"r")) as inp,\
            closing(openCompressed(tab_file,"w")) as out:
        w = None
        for line in inp:
            rec = annot_line_to_rec(line)
            if rec["recType"] != recType:
                continue
            if rec["recType"] == "annotPred":
                if int(rec["cntHits"]) < cntHits_min:
                    continue
            if int(rec["lenSamp"]) < lenSamp_min:
                continue
            if w is None:
                flds = sorted(rec.keys())
                w = csv.DictWriter(out, fieldnames=flds, restval='',dialect='excel-tab')
                w.writerow(dict([(fld,fld) for fld in flds]))
            w.writerow(rec)



def main():
    logging_config(detail="high")
    parser = ArghParser(description="Utilities for analysing output of annotation consensus classifier")
    parser.add_commands([
        as_tab
        ]
        )
    parser.dispatch()

if __name__ == "__main__":
    main()

