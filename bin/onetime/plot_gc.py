import sys
import numpy as np
from pyfasta import Fasta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    inpFasta = sys.argv[1]
except IndexError:
    print "Arguments: fasta_file"
    sys.exit(1)

# 100 bp window.
window = 100
fa = Fasta(inpFasta)

for seqid in fa.keys():
# get sequence as a numpy array with dtype='c'--char
    seq = np.array(fa[seqid], dtype='c')
    seq

    gcs = (seq == 'C') | (seq == 'G')
    gcs

# cast the booleans to ints.
    gcs = gcs.astype(np.uint8)
    gcs

    kern = np.ones(window)/ window
    kern

# same has boundary effects but output array is same length as seq
    gc_avg = np.convolve(gcs, kern, mode="same")
    _ = plt.figure(figsize=(8, 3))
    _ = plt.plot(gc_avg)
    _ = plt.ylim(0, 1)
    _ = plt.xlim(0, len(gc_avg))
    _ = plt.title("gc content for window:%i of :%s" % (window, seqid))
    plt.savefig('gc_%s_%i.png' % (seqid.split()[0],window))

