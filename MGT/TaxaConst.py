### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Various constants for taxonomy tree"""
#
## Names for various important taxonomy IDs
#

## @todo Phages vs other Viruses should be discriminated by divid: taxonomy/division.dmp

rootTaxid = 1
unclassRootTaxid = 12908
viralRootTaxid = 10239
virTaxid = viralRootTaxid
myovirTaxid = 10662 # Myoviridae - most of cyanophages are where
bacTaxid = 2
archTaxid = 2157
eukTaxid = 2759
micTaxids = (bacTaxid,archTaxid)
viroidsTaxid = 12884
cellTaxid = 131567 # cellular organisms - bact, arch & euk
micVirTaxids = (virTaxid,bacTaxid,archTaxid)
skingdomTaxids = (virTaxid,bacTaxid,archTaxid,eukTaxid,viroidsTaxid)

diatomsTaxid = 2836 # unicellular euk phytoplankton, major group of euk algae, divid plants
cyanobacTaxid = 1117 # Cyanobacteria phylum
shewanTaxid = 22 # Shewanella genus

viridiplantaeTaxid = 33090 #all green plants
chlorophytaTaxid = 3041 #green algae under viridiplantaeTaxid
streptophytaTaxid = 35493 #green plants - the only node other than chlorophytaTaxid under viridiplantaeTaxid, but
#it also includes some algae
embryophytaTaxid = 3193 #"land plants", "higher plants" under streptophytaTaxid, includes onion, rice etc
chordataTaxid = 7711 #euks, animals with chord, includes all vertebrates and some invertebrates

envBacTaxid = 48479 # Bacteria -> Environmental Samples node

viralTaxidLev2 = (\
        35237, #dsDNA
        35325, #dsRNA
        35268, #retroid
        29258, #ssDNA
        439488, #ssRNA
        )
phageTailedTaxid = 28883


## Main Linnaean ranks, modified with superkingdom in place of domain rank
linnMainRanks = ("species","genus","family","order","class","phylum","kingdom","superkingdom")
## We add a special linn rank for the root node, so that there is never a no_rank as a min
## rank in a given lineage. The root node rank must be modified when the tree is constructed.
## This list is used by TaxaLevels class.
linnMainRanksWithRoot = ("species","genus","family","order","class","phylum","kingdom","superkingdom","root")
## All Linnaenan ranks, modified as linnMainRanks
## @todo add sub- and super-ranks, such as subspecies and superfamily
linnRanks = linnMainRanks
viralRanksTemplate = linnMainRanks
noRank = "no_rank"
unclassRank = "unclassified"

## Override ranks for specific taxids when the tree is constructed
rankPatchesByTaxid = { 
        "root" : (rootTaxid,),
        "superkingdom" : (viralRootTaxid,) 
        }

## Some division IDs from NCBI division.dmp
dividEnv = 11

## We assert that NCBI taxonomy ID will never exceed this value.
##
## The reason for this to assertion is that we will assign our 
## own taxonomy IDs above that. Max ID as of 2010-10-15 is ~800K.
## Any estimates on the number of species out where are below 100M.
## @todo Change the taxonomy ID datatype to string and have a prefix
## denoting the namespace from which the ID originated (ncbi,mgt etc).
ncbiTaxidMax = 200000000

## Our own taxonomy IDs are generating starting from this
mgtTaxidFirst = ncbiTaxidMax + 100

