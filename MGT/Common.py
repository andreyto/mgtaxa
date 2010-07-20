### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from MGT.Types import *
from MGT.Util import *
from MGT.BatchRun import *
from MGT.Config import *
from MGT.UUID import *
from MGT.Debug import *

import numpy
import numpy.random as nrnd
n = numpy
nma = numpy.ma
from numpy import array
import random
from random import sample as sampleWOR #Random sampling w/o replacement
import math

def sampleBoundWOR(l,n):
    """Random sampling w/o replacement that clips sample size parameter by population size"""
    return sampleWOR(l,min(len(l),n))

def sampleRatioWOR(l,x,roundDir="round"):
    """Random sampling w/o replacement that tries to select a given ratio of the population.
    @param l population sequence
    @param x target ratio to select [0...1]
    @param roundDir rounding direction to get the integer number of samples to pick - one of 
    round, ceil or floor."""
    npop = len(l)
    ntarg = float(npop)*x
    if roundDir == "round":
        npick = int(round(ntarg))
    elif roundDir == "ceil":
        npick = int(math.ceil(ntarg))
    elif roundDir == "floor":
        npick = int(math.floor(ntarg))
    return sampleWOR(l,npick)

import os, sys, atexit, re, gzip
pjoin = os.path.join
from glob import iglob
import glob
import itertools
it = itertools
import operator
from types import StringTypes
from collections import defaultdict as defdict





