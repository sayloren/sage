"""

Script to produce a single interval file with the number of instances of overlap 
from a group of other bed files using bedtools intersect -c

Input: An interval file for your regions of interests, and a master txt file with 
the names of the other files to count intersects from

Output: An interval file with the original coordinates for your regions, 
with the number of interesections from your regions of interest AND an interval file
with the original coordiantes from your regions, with another column of the sum of 
the intersecting regions, ranked by most intersections

NOTE: will not take files with headers 

Version 2: headers added from file names

To Do: use bedtool annotate, also figure out how to convert to binary matrix, 
just in case

Wren Saylor
November 30 2016

"""

import argparse
import pybedtools as pbt
import pandas as pd


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

# 2 - turn into panda
def bedtoolToPanda(myRegions):
	matrixRegions = pd.read_table(myRegions.fn, header=None)
	return matrixRegions

# # 3 - convert to binary matrix (if were to do this, make new columns, then cut old)
# def binaryMatrix(matrixRegions):
# 	matrixRegions[''] = matrixRegions.iloc[:,4::].apply(lambda x: 1 if x >= 1 else 0)
# 	return 

# 3 - sum the columns
def pandaColSums(matrixRegions):
	matrixRegions['TotalCount'] = matrixRegions.iloc[:,4::].sum(axis=1) 
	pdRank = matrixRegions
	return pdRank

# 4 - sort by rank
def sortRank(pdRank):
	#pdSort = pdRank.sort(['TotalCount'],ascending=False)
	pdSort = pdRank.sort_values(['TotalCount'],ascending=False)
	return pdSort

# 5 - save file from panda
def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', index=False)


def main(args):

	myRegions = eachFileProcess(args.regions)
	# grab the file names from the master txt file
	aFiles = [line.strip() for line in args.file]
	myColumns = ['Chromosome','Start','Stop','ID']
	# for file in master txt file
	for fileName in aFiles:
		strFile = str(fileName)
		myColumns.append(strFile)
		btFeatures = eachFileProcess(fileName)
		myRegions = myRegions.intersect(btFeatures,c=True)
	matrixRegions = bedtoolToPanda(myRegions)
	matrixRegions.columns = myColumns
	pdRank = pandaColSums(matrixRegions)
	pdSort = sortRank(pdRank)
	savePanda(pdSort,"CountMatrix"+str(args.regions))


if __name__ == "__main__":
	args = get_args()
	main(args)