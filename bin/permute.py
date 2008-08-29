"""Read in a file, randomly permute the lines' order, write them out.
Loads an entire file into memory."""

import sys
import numpy

inp = open(sys.argv[1],'r')
print "Reading the input file %s" % sys.argv[1]
lines = inp.readlines()
inp.close()
n = len(lines)
ind = numpy.random.permutation(numpy.arange(n))
out = open(sys.argv[2],'w')
print "Writing the output file %s" % sys.argv[2]
for i in ind: out.write(lines[i])
out.close()
print "Done, wrote %s lines" % n

