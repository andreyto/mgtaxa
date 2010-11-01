### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Service script that inserts the copyright reference into every source file in the directory tree"""
import os
pjoin = os.path.join
from copy import copy
import sys

topSrcDir = os.environ["MGT_HOME"]
insLines = open(pjoin(topSrcDir,"etc/source_header"),'r').readlines()
insLines.append("\n")
# make language-specific comments in the lines to be inserted
# original lines already start with python/shell comments
insLinesPy = insLines
insLinesC = [ "//"+l for l in insLines ]

iFile = 0
for (d,ds,fs) in os.walk(topSrcDir):
    ds_cp = copy(ds)
    for x in ds_cp:
        if x.startswith("."):
            ds.remove(x)
        elif d == topSrcDir and x in ("build","external","test_external","doc"):
            ds.remove(x)
    for f in fs:
        if f.endswith(".py"):
            ins = insLinesPy
        elif sum((f.endswith(ext) for ext in (".c",".cpp",".h",".hpp"))):
            ins = insLinesC
        elif pjoin(d,f).startswith(pjoin(topSrcDir,"bin")+'/'):
            ins = insLinesPy
        else:
            ins = None
        if ins is not None:
            srcFile = pjoin(d,f)
            inp = open(srcFile,'r')
            lines = inp.readlines()
            inp.close()
            iFile += 1
            skip = False
            for l in lines[:10]:
                if "COPYING" in l:
                    print "%s done before, skipping" % srcFile
                    skip = True
            if not skip:
                print srcFile, ins[0]
                if len(lines) > 0 and lines[0].startswith("#!"):
                    lines = lines[0:1] + ins + lines[1:]
                else:
                    lines = ins + lines
                out = open(srcFile,'w')
                for l in lines:
                    out.write(l)
                out.close()


