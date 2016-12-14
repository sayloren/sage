"""

Script for converting density data from UCSC (was wig, used wig2bed, now bedfile where 
fourth column is density) to the bin-ed density needed to run in Ruth's Density script

Input: A master text file with the file names of the interval files you wish to convert
and the bin file returned from Ruth's Density script with the windows to use

Output: Density files with proper windows binned that can be piped right into Ruth's Matrix
maker

Version 2: added sort

Wren Saylor, adapted from Ruth's Density Script
December 13 2016

"""

import argparse
import pybedtools as pbt
import pandas as pd


def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a'
	'list of paths to the feature files you want to process, '
	'separated by newlines')
	parser.add_argument("-b", "--bin", type=str, default="Bins_for_density_0based.bed")

	return parser.parse_args()

# 1 - read in files
def getFileNames(args):
	arFilenames = [line.strip() for line in args.file]
	return arFilenames

# 2 - get features of file
def getFeatures(strFileName):
	preSortbtFeatures = pbt.BedTool(strFileName)
	btFeatures = preSortbtFeatures.sort()
	return btFeatures
	
# 3 - left outer join intersect 
def intersectWindow(btWindows, btFeatures):
	btIntersect = btWindows.intersect(btFeatures,loj=True)
	return btIntersect

# 4 - return windows with attached densities
def cutToFraction(btIntersect):
	btCutToFraction = btIntersect.cut([0,1,2,6])
	return btCutToFraction

# 5 - perform the 3 and 4
def getDensity(btWindows, btFeatures):
	btIntersect = intersectWindow(btWindows, btFeatures)
	btCutToFraction = cutToFraction(btIntersect)
	return btCutToFraction

# 6 - intersect, then merge the features in the window
def mergeWindow(btCutToFraction):
	btMerge = btCutToFraction.merge(c=4,o="mean")
	return btMerge

# 7 - save to bedtool
def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

# 8 - string commands together
def forEachFeatureFile(strFileName, btWindows):
	btFeatures = getFeatures(strFileName)
	btCutToFraction = getDensity(btWindows, btFeatures)
	btMerge = mergeWindow(btCutToFraction)
	saveBedTool(btMerge, str('Density_for_matrix_{0}.bed'.format(strFileName)))
	return btMerge

def main():
	args = get_args()
	btWindows = getFeatures(args.bin)
	arFilenames = getFileNames(args)
	for filename in arFilenames:
		btMerge = forEachFeatureFile(filename, btWindows) #, btCutToFraction

if __name__ == "__main__":
	main()

