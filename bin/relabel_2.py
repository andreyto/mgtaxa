### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


inp = open("test.svm",'r')
out = open("test.svm.1",'w')

for line in inp:
    w = line.split(None,1)
    if w[0] == '1':
        l = '15'
    elif w[0] == '2':
        l = '16'
    else:
        raise ValueError(w[0])
    out.write("%s %s" % (l,w[1]))
out.close()

