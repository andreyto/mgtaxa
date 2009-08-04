### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


tr = open("train.svm",'r').readlines()
ts = open("test.svm",'r').readlines()
#trv = [ x.split(None,1)[1] for x in tr ]
#tsv = [ x.split(None,1)[1] for x in ts ]
trv = tr
tsv = ts
trs = set(trv)
tss = set(tsv)
tss_i_trs = tss.intersection(trs)
print "len(ts) = %s len(tss) = %s len(tr) = %s len(trs) = %s len(tss_i_trs) = %s" % \
        (len(ts),len(tss),len(tr),len(trs),len(tss_i_trs))
tss_d_trs = tss.difference(trs)
ts_out = open('test.clean.svm','w')
ts_out.write(''.join(sorted(list(tss_d_trs))))
ts_out.close()

