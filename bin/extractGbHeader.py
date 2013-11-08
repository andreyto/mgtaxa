### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Extract everything from GenBank stream but the ORIGIN sequence.
Reads standard input, writes to standard ouput"""

import sys
from MGT.Common import *

skip = False
nSeqLines = 0

try:
    inp = openCompressed(sys.argv[1],"r")
except IndexError:
    inp = sys.stdin

try:
    out = openCompressed(sys.argv[2],"w")
except IndexError:
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

