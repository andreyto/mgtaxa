from MGT.Common import *
from MGT.Logging import *
from MGT.FastaIO import *

from argh import ArghParser,arg

log = logging.getLogger(os.path.basename(sys.argv[0]))

def rec_filter(inp_file,out_file,len_min=0,len_max=10**18,out_line_len=80,id_list=None):
    """Filter records from a FASTA file"""
    if id_list is not None:
        with open(id_list,"r") as inp_id_list:
            id_list = set(( line.strip() for line in inp_id_list))
    with closing(FastaReader(inp_file)) as inp,\
            closing(FastaWriter(out_file)) as out:
        for rec in inp.records():
            hdr = rec.header()
            seq = rec.sequence()
            if len(seq) < len_min or len(seq) > len_max:
                continue
            if id_list is not None and rec.getId() not in id_list:
                continue
            out.record(hdr,seq)


def main():
    logging_config(detail="high")
    parser = ArghParser(description="Utilities for manipulating FASTA files")
    parser.add_commands([
        rec_filter
        ]
        )
    parser.dispatch()

if __name__ == "__main__":
    main()

