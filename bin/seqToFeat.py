### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Read a file with sequence features and create a file with sparse real features based on word distance histogram."""

from MGT.App import *
from MGT.SeqFeaturesApp import SeqFeaturesApp

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(SeqFeaturesApp)

