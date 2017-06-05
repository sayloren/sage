"""

Script for converting density data from UCSC (was wig, used wig2bed, now bedfile where 
fourth column is density) to the bin-ed density needed to run in Ruth's Density script

Input: A master text file with the file names of the interval files you wish to convert
and the bin file returned from Ruth's Density script with the windows to use

Output: Density files with proper windows binned that can be piped right into Ruth's Matrix
maker

To Do: The merge might be over merging the lines - compressing instead of bin-ing

Version 2: added sort

Version 3: streamlined function

Version 4: used alternative to bedtools merge for combining on density, to avoid over merging

Version 5: used map instead of intersect + merge to get the density over the input bins,
also allows easier use of mode argument

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
	parser.add_argument("-m", "--mode", type=str, default="mean", help='the method of'
	'selecting how to get the density value for the converted file '
	'look at http://bedtools.readthedocs.io/en/latest/content/tools/map.html for options')
	return parser.parse_args()

# 1 - read in files
def getFileNames(args):
	arFilenames = [line.strip() for line in args.file]
	return arFilenames

# 2 - get features of file
def getFeatures(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	#btFeatures = preSortbtFeatures.sort()
	return btFeatures

# 3i - left outer join intersect - depreciated; will use map instead of intersect
def intersectWindow(btWindows, btFeatures):
	btIntersect = btWindows.intersect(btFeatures,loj=True)
	return btIntersect

# 3 - map density of your input feature files to the bins
def mapWindow(btWindows,btFeatures,mode):
	btMap = btWindows.map(btFeatures,c=4,o=mode)

# 4 - grab just columns for window coordinates and densities - depreciated; used in conjunction with intersect
def cutToFraction(btIntersect):
	btCutToFraction = btIntersect.cut([0,1,2,6])
	return btCutToFraction

# 5 - bedtool to panda for merge to be implemented correctly - depreciated; used in conjunction with intersect
def bedtoolToPanda(btobject):
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

# 5i - merge the density features in the windows by mean - depreciated; will merge more extensively than desired
def mergeWindow(btCutToFraction,mode):
	btMerge = btCutToFraction.merge(c=4,o=mode)
	return btMerge

# 6 - merge density features - depreciated; used in conjunction with intersect
def mergeWindow(pdFraction): #,mode
	pdFraction.columns = ["chr","start","stop","den"]
	pdFraction.loc[pdFraction.den == '.'] = 0
	pdFraction[["den"]] = pdFraction[["den"]].astype(float)
	pdGroup = pdFraction.groupby(["chr","start","stop"],as_index=False)["den"].mean()
	return pdGroup

# 7i - save to bedtool - depreciated; no longer output from bedtool
def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

# 7 - save file from panda
def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', header=False, index=False)

def main():
	args = get_args()
	btWindows = getFeatures(args.bin)
	arFilenames = getFileNames(args)
	for filename in arFilenames:
		btFeatures = getFeatures(filename)
		btMap = mapWindow(btWindows,btFeatures,args.mode)
# 		btIntersect = intersectWindow(btWindows,btFeatures)
# 		btCutToFraction = cutToFraction(btMap)
# 		pdFraction = bedtoolToPanda(btCutToFraction)
# 		btMerge = mergeWindow(pdFraction)#,args.mode
		savePanda(btMap, str('Density_for_matrix_{0}.bed'.format(filename)))

if __name__ == "__main__":
	main()

