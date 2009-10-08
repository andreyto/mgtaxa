### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from MGT.Util import *

inpFile = pjoin(options.testDataDir,"iofilt.test.gz")
inpFile = "/usr/local/projects/GOSII/syooseph/MF150/all_Moore_data/CP000878.1.gbk.gz"
inp = ioFilter(openCompressed(inpFile,'r'),code="lambda x: x.replace('&gt;','')",mode="line")

ok = False
for line in inp:
    print line,
    if "complement(1356723..1357840)" in line:
        ok = True
assert ok

