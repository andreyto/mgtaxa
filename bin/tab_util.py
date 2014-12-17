from MGT.Common import *
from MGT.Logging import *

import itertools, csv, collections
from collections import OrderedDict
from argh import ArghParser,arg

log = logging.getLogger(os.path.basename(sys.argv[0]))

@arg("inp-file",help="Input file")
@arg("out-file",help="""Output file. If the
    first non-skipped result is a mapping type, output will be
    formatted with csv.DictWriter and start with a header, otherwise, 
    it must be a sequence type and output will be formatted with 
    csv.writer.""")
@arg("--expr",help="""Python expression that will be executed with Python
    `eval` statement and given a locals() envronment containing
    a variable called `x` set to the current record parsed either
    with csv.DictReader or csv.reader. Expression evaluation will have access to OrderedDict
    definition if the order of columns is important.""")
@arg("--where",help="""Expression with an argument `x` where `x` is a result of --expr;
    when evaluating to False is a logical scope, the record will be skipped from
    the output""")
@arg("--fields",help="Names or positions of fields to select from `expr` output")
@arg("--inp-header",help="If True, assume the first input row contains column names")
@arg("--no-out-header",help="If True, output column names if --inp-header was also True")

def filter_rows(inp_file,
        out_file,
        expr='x',
        where='True',
        fields=None,
        inp_header=True,
        no_out_header=False):
    """Transform/filter each row of the input tabular file"""
    out_is_dict = None
    fields_filt = lambda y: y
    with closing(openCompressed(inp_file,"r")) as inp,\
            closing(openCompressed(out_file,"w")) as out:
        if inp_header:
            r = csv.DictReader(inp,dialect='excel-tab')
        else:
            r = csv.reader(inp,dialect='excel-tab')
        w = None
        for inp_rec in r:
            out_rec = eval(expr,{},dict(x=inp_rec))
            where_res = eval(where,{},dict(x=out_rec))
            if not where_res:
                continue
            if w is None:
                if isinstance(out_rec, collections.Mapping):
                    out_is_dict = True
                else:
                    out_is_dict = False
                if fields is not None:
                    fields = fields.strip().split(",")
                    if out_is_dict:
                        fields_filt = lambda y: dict(((f,y[f]) for f in fields))
                    else:
                        fields = [ int(f) for f in fields ]
                        fields_filt = lambda y: [ y[f] for f in fields ]
                out_rec = fields_filt(out_rec)
                if out_is_dict:
                    flds = out_rec.keys()
                    w = csv.DictWriter(out, fieldnames=flds, restval='',dialect='excel-tab')
                    if not no_out_header:
                        w.writerow(dict([(fld,fld) for fld in flds]))
                else:
                    w = csv.writer(out,dialect='excel-tab')
            else:
                out_rec = fields_filt(out_rec)
            w.writerow(out_rec)



def main():
    logging_config(detail="high")
    parser = ArghParser(description="Utilities for working with tabular files")
    parser.add_commands([
        filter_rows
        ]
        )
    parser.dispatch()

if __name__ == "__main__":
    main()

