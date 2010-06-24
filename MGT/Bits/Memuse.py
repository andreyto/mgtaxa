"""
Get memory usage for the process.
Adapted from http://code.activestate.com/recipes/286222/
"""
# A.T.
# Added string formatted output.
# Notes on the alternatives:
# resource.getrusage() returns zeros because most fields are
# not implemented in 2.6 Linux kernel.
# Other alternatives would be to run Linux utilities like 
# pmem (that gives the total on virtual memory on the last line of output)

import os

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}
_scale_inv = ( (1024.*1024.,"MB"),(1024.,"KB") )

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]

def toScale(x):
    for sc in _scale_inv:
        y = x/sc[0]
        if y >= 1:
            return "%.3f%s" % (y,sc[1])
    return "%.3f%s" % (y,"B")

def memory(since=0.0,asStr=True):
    '''Return memory usage in bytes or as formatted string.
    '''
    b = _VmB('VmSize:') - since
    if asStr:
        return "VirtMem: " + toScale(b)
        


def resident(since=0.0,asStr=True):
    '''Return resident memory usage in bytes.
    '''
    b = _VmB('VmRSS:') - since
    if asStr:
        return "ResMem: " + toScale(b)


def stacksize(since=0.0,asStr=True):
    '''Return stack size in bytes.
    '''
    b = _VmB('VmStk:') - since
    if asStr:
        return "StackMem: " + toScale(b)

