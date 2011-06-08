"""Get per-contig read counts from 454ReadStatus.txt produced by Newbler assembler.
We count 5' and 3' ends as 1/2 each to handle split reads."""

import sys
from collections import defaultdict as defdict

cnt = defdict(float)

# skip header
sys.stdin.next()
# count contigs, including cases when read spans two contigs
# (we assign 1/2 to each - what else to do?)
# FRDNTY201A47E7  Assembled       contig208724    402     -       contig134302    1772    +
for line in sys.stdin:
    words = line.strip().split()
    if len(words) >= 3:
        cnt[words[2]] += 0.5
        cnt[words[5]] += 0.5

out = sys.stdout
for (contig,count) in sorted(cnt.iteritems()):
    out.write("%s\t%i\n" % (contig,round(count)))

