"""

Script to produce a single interval file with the number of instances of overlap 
from a group of other bed files using bedtools intersect -c

Input: An interval file for your regions of interests, and a master txt file with 
the names of the other files to count intersects from

Output: An interval file with the original coordinates for your regions, 
with the number of interesections from your regions of interest

NOTE: will not take files with headers - and will not annotate in headers with printing
order of files in master txt file will be order of columns of intersections to your regions
This could be done by converting to a pandas data frame: To Do v.2

Wren Saylor
November 30 2016

"""

import argparse
import pybedtools as pbt


def get_args():
	parser = argparse.ArgumentParser(description='Script to produce a single interval '
		'file with the number of instances of overlap from a group of other bed files '
		'using bedtools intersect -c')
	parser.add_argument("-r", "--regions", type=str, required=True, help= "Interval file for your coordiantes") #argparse.FileType('rU')
	parser.add_argument("file", type=argparse.FileType('rU'), help= "Regions to count overlaps from")
	return parser.parse_args()

# 1 - get the features for each file
def eachFileProcess(fileName):
	btFeatures = pbt.BedTool(fileName)
	return btFeatures

# 2 - save the file as slop_filename.bed
def saveBedTool(xxx, strFilename):
	xxx.saveas(strFilename)


def main(args):

	myRegions = eachFileProcess(args.regions)

	# grab the file names from the master txt file
	aFiles = [line.strip() for line in args.file]
	# for file in master txt file
	for fileName in aFiles:
		btFeatures = eachFileProcess(fileName)
		myRegions = myRegions.intersect(btFeatures,c=True)

	saveBedTool(myRegions,"CompiledCounts"+str(args.regions))


if __name__ == "__main__":
	args = get_args()
	main(args)