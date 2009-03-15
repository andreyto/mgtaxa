"""Various constants for taxonomy tree"""
#
## Names for various important taxonomy IDs
#

## @todo Phages vs Viruses should be discriminated by divid: taxonomy/division.dmp

unclassRootTaxid = 12908
viralRootTaxid = 10239
virTaxid = viralRootTaxid
myovirTaxid = 10662 # Myoviridae - most of cyanophages are where
bacTaxid = 2
archTaxid = 2157
eukTaxid = 2759
micTaxids = (bacTaxid,archTaxid)

diatomsTaxid = 2836 # unicellular euk phytoplankton, major group of euk algae, divid plants
cyanobacTaxid = 1117 # Cyanobacteria phylum
shewanTaxid = 22 # Shewanella genus

viridiplantaeTaxid = 33090 #all green plants
chlorophytaTaxid = 3041 #green algae under viridiplantaeTaxid
streptophytaTaxid = 35493 #green plants - the only node other than chlorophytaTaxid under viridiplantaeTaxid, but
#it also includes some algae
embryophytaTaxid = 3193 #"land plants", "higher plants" under streptophytaTaxid, includes onion, rice etc
chordataTaxid = 7711 #euks, animals with chord, includes all vertebrates and some invertebrates

viralTaxidLev2 = (\
        35237, #dsDNA
        35325, #dsRNA
        35268, #retroid
        29258, #ssDNA
        439488, #ssRNA
        )
phageTailedTaxid = 28883
viroidsTaxid = 12884
cellTaxid = 131567 # cellular organisms - bact, arch & euk
micVirTaxids = (virTaxid,bacTaxid,archTaxid)
skingdomTaxids = (virTaxid,bacTaxid,archTaxid,eukTaxid,viroidsTaxid)


## Main Linnaean ranks, modified with superkingdom in place of domain rank
linnMainRanks = ("species","genus","family","order","class","phylum","kingdom","superkingdom")
viralRanksTemplate = linnMainRanks
noRank = "no_rank"
unclassRank = "unclassified"

