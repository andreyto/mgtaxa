import sys
import itertools

if len(sys.argv) < 2:
    print "Usage: {0} min_contig_length < annot.csv > annot.filtered.csv"
    sys.exit(1)

minLenSamp = int(sys.argv[1])

for line in sys.stdin:
    words = line.split("\t")
    map = dict(itertools.izip(words[::2],words[1::2]))
    if int(map.get("lenSamp",0)) >= minLenSamp:
        sys.stdout.write(line)

