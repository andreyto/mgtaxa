### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Extract everything from GenBank stream but the ORIGIN sequence.
Reads standard input, writes to standard ouput"""

import sys

skip = False
nSeqLines = 0
inp = sys.stdin
out = sys.stdout
while True:
    line = inp.readline()
    if not line:
        break
    if line.startswith('ORIGIN'):
        out.write(line)
        line = inp.readline()
        out.write(line)
        if not line.startswith('//'):
            while True:
                line = inp.readline()
                if line.startswith('//'):
                    out.write(line)
                    break
    else:
        out.write(line)

