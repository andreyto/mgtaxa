from MGT.Common import *


def test_runsOfOnesArray():
    data = [ array([0,0,0,1,1,1,0,0,0,1,1,1,1,0,0]),
        array([1,0,0,1,1,1,0,0,0,1,1,1,1,0,0]),
        array([1,0,0,1,1,1,0,0,0,1,1,1,1,0,1]),
        array([0,0,0,1,1,1,0,0,0,1,1,1,1,0,1]),
        array([0,0,0]),
        array([1,1,1]),
        array([1]),
        array([0]), 
        array([]) ]

    for b in data:
        print "Input: ", b
        runs = runsOfOnesArray(b)
        print "Output: ", runs
        print "Input[Output]:"
        s = 0
        for r in runs:
            v = b[r[0]:r[1]]
            print " "*4,v
            assert n.alltrue(v)
            s += v.sum()
        assert b.sum() == s
        print "*****************\n"

test_runsOfOnesArray()

