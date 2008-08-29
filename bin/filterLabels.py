"""Filter SVMLight file by a label set."""

import sys

iArgv = 1
indexOutFile = sys.argv[iArgv]
iArgv+=1

iArgv = 1
labels = set([ int(x) for x in sys.argv[iArgv:] ])

for line in sys.stdin:
    lab, feature = line.split(None,1)
    if int(lab) in labels:
        sys.stdout.write(line)

