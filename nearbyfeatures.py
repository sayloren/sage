"""

Script to look at how close features are between two bed files, or a UCE file and a list of bedfiles

Wren Saylor
Jan 2018

To Do:
-Stats per bed file features (number, size, genome coverage)
-Stats per overlaps (min/max, average, std) 
-Where do overlaps occur

Copyright 2018 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""
import argparse
import pandas as pd
import pybedtools as pbt

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the tertiary features to query")# Genes
	parser.add_argument("-g","--genomefile",type=str,help="genome file",default='hg19.genome')
	return parser.parse_args()

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get stats
def descriptivestats(btfeature,genomefile):
# 	coverage = btfeature.genome_coverage(genome=genomefile) # genome coverage, not working yet
	totalcount = btfeature.count() # total number of features in bedfile
	
	
	
	# distance to other features in group
	# size (ave, std, range) - as pd
	# get relevant info into pd
	return pdstats



def main():
	args = get_args()
	
	# 1) read in files from args
	primaryfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tertiaryfiles = [line.strip() for line in args.tertiaryfeatures]
	
	genomefile = args.genomefile
	
	allfiles = secondaryfiles + tertiaryfiles
	allfiles.insert(0,primaryfile)
	
	print 'running script for {0}'.format(allfiles)
	
	# 2) descriptive stats per feature
	for file in allfiles:
		btfeature = get_bedtools_features(file) # get bed file features
		btstats = descriptivestats(btfeature,genomefile) # get some stats about each file
		# collect stats into df
		# save stats file
	
	
	# 3) get overlaps - secondary; intersect
	# 4) get nearest - tertiary; closest



if __name__ == "__main__":
	main()