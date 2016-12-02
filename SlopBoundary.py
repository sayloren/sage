"""

Script to produce interval files for varying slopped edges, 
sorting and mergeing the interval files after slopping, and before writing to file

Input: A file with the file names of 0_based interval files

OutputL: A new file for each file in master txt file x slop size combination, labeled
in format slopsize#_filename.bed

Wren Saylor 
November 30 2016
 
"""

import argparse
import pybedtools as pbt


def get_args():
	parser = argparse.ArgumentParser(description='Script to produce interval files '
		'for varying slopped edges, sorting and mergeing the interval files after slopping, '
		'and before writing to file')
	parser.add_argument("file", type=argparse.FileType('rU'))
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
	parser.add_argument("-s", "--slop", type=int, help='Size of desired slopped '
		'edges separated by spaces', nargs='+')
	return parser.parse_args()

# 1 - get the features for each file
def eachFileProcess(fileName):
	btFeatures = pbt.BedTool(fileName)
	return btFeatures

# 2 - slop the edges
def getSlopped(btFeatures,aSize,aGenome):
	btSlopped = btFeatures.slop(b=aSize, g=aGenome)
	print btSlopped
	return btSlopped

# 3 - sort and merge
def sortAndMerge(btSlopped):
	btSorted = btSlopped.sort()
	btMerged = btSorted.merge()
	return btMerged


# 4 - save the file as slop_filename.bed
def saveBedTool(btMerged, strFilename):
	btMerged.saveas(strFilename)


def main(args):

	# grab the file names from the master txt file
	aFiles = [line.strip() for line in args.file]
	
	# put the genome file into a string
	aGenome = args.genome
	
	# for file in master txt file
	for fileName in aFiles:
	
		# for slop size in input slop sizes
		for aSize in args.slop:
			btFeatures = eachFileProcess(fileName)
			btSlopped = getSlopped(btFeatures,aSize,aGenome)
			btMerged = sortAndMerge(btSlopped)
			saveBedTool(btMerged,str(aSize) + "_" + str(fileName))
			

if __name__ == "__main__":
	args = get_args()
	main(args)
