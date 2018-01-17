"""

Script to look at how close features are between two bed files, or a UCE file and a list of bedfiles

Wren Saylor
Jan 2018

To Do:
-Stats per bed file features (number, size, genome coverage)
-Stats per overlaps (min/max, average, std) 
-Where do overlaps occur

"""
import argparse
import pandas as pd
import pybedtools as pbt

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the tertiary features to query")# Genes
	parset.add_argument("-g","--genome",type=argparse.FileType('rU'),help="genome file",default='hg19.genome')
	return parser.parse_args()

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get stats
def descriptivestats(btfeature,genomefile):
	genomecoverage = btfeature.genome_coverage(g=genomefile) # genome coverage
	totalcount = btfeature.count() # total number of features in bedfile
	# size (ave, std, range)
	# distance to other features in group
	# get relevant info into pd
	return pdstats



def main():
	args = get_args()
	
	# 1) read in
	primaryfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tertiaryfiles = [line.strip() for line in args.tertiaryfeatures]
	
	genome = args.genome
	
	allfiles = primaryfile + secondaryfiles + tertiaryfiles
	
	# 2) descriptive stats per feature
	for file in allfiles:
		btfeature = get_bedtools_features(file)
		btstats = descriptivestats(btfeature,genome)
		# collect stats into df
		# save stats file
	
	
	# 3) get overlaps - secondary; intersect
	# 4) get nearest - tertiary; closest



if __name__ == "__main__":
	main()