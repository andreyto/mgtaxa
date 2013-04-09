"""Compute coverage stats for Celera Assembler contigs"""
from MGT.Common import *

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from joblib import Memory
makedir("cache.tmp")
memory = Memory(cachedir="cache.tmp", verbose=0)


@memory.cache
def ctg_cov(frg_ctg_files,ctg_len_files,label):

    ctg_frg_bases = defdict(int)
    ctg_frgs = defdict(int)

    for frg_ctg in frg_ctg_files:
        inp = open(frg_ctg,'r')
        for line in inp:
            words = line.split()
            ctg_frg_bases[words[1]] += abs(int(words[3])-int(words[2]))
            ctg_frgs[words[1]] += 1
        inp.close()

    ctg_frg_bases = dict(ctg_frg_bases)
    ctg_frgs = dict(ctg_frgs)

    ctgs = n.zeros(len(ctg_frg_bases),
            dtype=[("id",idDtype),
                ("frg_bases",int),
                ("frgs",int),
                ("ctg_len",int),
                ("cov",float),
                ("label",int)])

    i_ctg = 0
    for (ctg_len,lab) in it.izip(ctg_len_files,label):
        inp = open(ctg_len,'r')
        for line in inp:
            try:
                words = line.split()
                id_ctg = words[0]
                ctg = ctgs[i_ctg]
                ctg["id"] = id_ctg
                ctg["frg_bases"] = ctg_frg_bases[id_ctg]
                ctg["frgs"] = ctg_frgs[id_ctg]
                ctg["ctg_len"] = int(words[1])
                ctg["cov"] = float(ctg["frg_bases"])/ctg["ctg_len"]
                ctg["label"] = lab
                i_ctg += 1
            except KeyError,key:
                print "No ctg_frg record for id %s" % (key,)
        inp.close()
    assert i_ctg == len(ctgs)

    return ctgs

labCtg = 1
labDeg = 2

#/usr/local/archive/projects/GOS/OMZ/Assembly/CA/9-terminator

term_dir = sys.argv[1]
asm_pref = sys.argv[2]
qc_TotalUsableReads=int(sys.argv[3])

frg_ctg = pjoin(term_dir,asm_pref+".posmap.frgctg")
frg_deg = pjoin(term_dir,asm_pref+".posmap.frgdeg")
ctg_len = pjoin(term_dir,asm_pref+".posmap.ctglen")
deg_len = pjoin(term_dir,asm_pref+".posmap.deglen")

ctgs = ctg_cov(frg_ctg_files=(frg_ctg,frg_deg),
        ctg_len_files=(ctg_len,deg_len),
        label=(labCtg,labDeg))

#plt.hist(ctgs["cov"],50)
max_len = ctgs["ctg_len"].max()
plt.hist(ctgs["cov"],100,range=(1,20),
        weights=ctgs["ctg_len"]/float(max_len),
        normed=1)
plt.xlabel("Coverage")
plt.ylabel("Bases")
plt.savefig("hist_bases_cov.png")
plt.close()

#This is a lot of points - takes a while
plt.scatter(ctgs["cov"],ctgs["ctg_len"])
plt.xlabel("Coverage")
plt.ylabel("Contig length")
plt.savefig("scat_cov_len.png")
plt.close()

ctgs_sel = ctgs[ctgs["cov"]<100]
plt.scatter(ctgs_sel["cov"],ctgs_sel["ctg_len"])
plt.xlabel("Coverage")
plt.ylabel("Contig length")
plt.savefig("scat_cov_100_len.png")
plt.close()

ctgs_sel = ctgs[n.logical_and(ctgs["cov"]<100,ctgs["ctg_len"]>5000)]
plt.scatter(ctgs_sel["cov"],ctgs_sel["ctg_len"])
plt.xlabel("Coverage")
plt.ylabel("Contig length")
plt.savefig("scat_cov_100_len_5000.png")
plt.close()

ctgs_sel = ctgs[n.logical_and(ctgs["cov"]<90,ctgs["cov"]>55)]
print "Avg frg size in selection: ", 
    float(ctgs_sel["frg_bases"].sum())/ctgs_sel["frgs"].sum()

ctgs_sel = ctgs[n.logical_and(ctgs["cov"]<18,ctgs["cov"]>0)]
print "Avg frg size in selection: ", 
    float(ctgs_sel["frg_bases"].sum())/ctgs_sel["frgs"].sum()

ctgs_sel = ctgs[n.logical_and(ctgs["cov"]<18,ctgs["ctg_len"]>10000)]
print "Avg frg size in selection: ", 
    float(ctgs_sel["frg_bases"].sum())/ctgs_sel["frgs"].sum()

cnt_len, bins = n.histogram(ctgs["cov"],100,range=(1,200),weights=ctgs["ctg_len"])
cnt_len_sq, bins = n.histogram(ctgs["cov"],100,range=(1,200),weights=ctgs["ctg_len"]**2)
cnt, bins = n.histogram(ctgs["cov"],bins)
avg_ctg_len = cnt_len.astype(float)/cnt
std_ctg_len = (cnt_len_sq.astype(float)/cnt - avg_ctg_len**2)**0.5
plt.bar(bins[:-1],avg_ctg_len,yerr=std_ctg_len)
plt.xlabel("Coverage")
plt.ylabel("Avg Contig length and Std")
plt.savefig("hist_cov_len.png")
plt.close()

#CA coverage stats (A-statistics)
#http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=Celera_Assembler_Theory
#http://www.sciencemag.org/content/287/5461/2196.full#ref-24
#"""
#Suppose there are F fragments in a database and the genome size is estimated as G. 
#For a unitig with k fragments and distance p between the start of its first fragment
#and the start of its last fragment, the probability of seeing the k - 1 start points
#in the interval of length p, given the unitig is not oversampled, is 
#[(pF/G)^k/k!]exp(-pF/G). If the unitig was the result of collapsing two repeats, then
#the probability is [(2pF/G)^k/k!]exp(-2pF/G). The log of the ratio of these two 
#probabilities, (log e)pF/G - (log 2) k, is our A statistic.
#"""
#The math for the expected length of contigs etc is in
#Arratia R, Lander E, Tavare S, Waterman M (1991) Genomic mapping by anchoring 
#random clones: A mathematical analysis. Genomics 11: 806-827
#http://www.cmb.usc.edu/papers/msw_papers/msw-098.pdf
#Math for metagenomics:
#Stanhope SA (2010) Occupancy Modeling, Maximum Contig Size Probabilities and 
#Designing Metagenomics Experiments. PLoS ONE 5(7): e11652. doi:10.1371/journal.pone.0011652
#CA reports computed genome size here: 5-consensus-coverage-stat/<asm_pref>.stats

#We want to make sure that even fairly high coverage bugs are not
#marked as degenerates. We use contigs that are not too short.
ctgs_sel =  ctgs[n.logical_and(ctgs["ctg_len"]>1000,ctgs["cov"]>50)]
print "Est genome length:", float(qc_TotalUsableReads)/ctgs_sel["frgs"].sum()*ctgs_sel["ctg_len"].sum()
