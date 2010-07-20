"""Datatypes and datatype sizes common among multiple modules"""

## NCBI taxid is an integer of max this length
taxidBytes = 4 

## Numpy dtype string for NCBI taxid
taxidDtype = "i%i" % taxidBytes

##@todo convert all code the uses hard-wired taxid size into using the vars from this module

