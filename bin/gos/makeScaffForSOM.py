"""Collects scaffolds from various GOS assemblies and combine them with chunks of NCBI RefSeq for SOM building.
The SOM should be then trained by whatever external method.
The resulting SOM should be loaded by a loader MrMpiSOM that should create a SOMModel object,
which is then used to map the input vectors, and then graphed by SOMGraph module.
"""
# we need to label separately these things (>=5Kb):
# NCBI bacteria Moore genomes
# NCBI bacteria non-Moore genomes
# NCBI phages
# NCBI non-phage viruses
# GOS phages from viral fraction
# GOS phages reliably classified as such by APIS from large fraction
# GOS scaffolds reliably classified as non-phage by APIS
