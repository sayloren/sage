# README #

### This is a handful of unrelated sripts that have been useful for UCE related things in the Wu Lab ###

IntersectFlagC.py 
Inputs: A file of regions of interests, a master.txt file of things to count overlaps for
Outputs: A file of regions of interest with the number of times the interval files in master.txt intersect

SlopBoundary.py
Inputs: A master.txt file with the interval filenames you want slopped and the size of the slops you want
Outputs: A new file for each file sloped and each size slop you requested

ScaleUseCounts.py
Inputs: A master.txt file with interval files, with a fourth column for rank
Outputs: A new interval file per input file from master.txt list, scaled to the range used in the UCSC table browser so that it can be viewed with those elements with a higher score receiving a darker shade of labeling

Unique.py
Inputs: The newly created UCEs and the comprehensive list of all old UCEs, both in fasta format
Outputs: The new and old UCE being compared, and whether or not they are the same


DensityConverter.py
Inputs: A master.txt file with the interval files, with the fourth density column, and the bin file from Ruth's density script
Outputs: A file for each input file, where the density is bined into the proper bins from Ruth's script

###
