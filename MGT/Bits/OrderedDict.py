"""Drop-in replacement for collections.OrderedDict for python < 2.7.
This will import collections.OrderedDict if available.
From http://code.activestate.com/recipes/576693/
"""

__all__ = [ "OrderedDict" ] 

try:
    from collections import OrderedDict
except ImportError: # version < python2.7
    from _OrderedDictImpl import *
 
