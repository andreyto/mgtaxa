### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Redirects to SeqFeaturesApp scripting interface"""

from MGT.App import *
from MGT.SeqFeaturesApp import SeqFeaturesApp

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(SeqFeaturesApp)

