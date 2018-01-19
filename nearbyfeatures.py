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

# convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	pdfeature = pd.read_table(btfeature.fn,header=None)
	return pdfeature

# label chr start and stop, make size
def label_coordinate_columns(pdfeature):
	pdfeature['size'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	return pdfeature

# get stats
def descriptive_stats(pdfeature,column):
	pdstat = pdfeature[column].describe()
	return pdstat

# save panda
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	intersect = afile.intersect(bfile,c=True)
	return intersect

def count_number_with_zero_overlaps(df,column):
	print len(df[(df[column]==0)])

# move elements without any overlaps
def remove_rows_with_no_overlaps(overlaps,column):
	pdfeature = overlaps[overlaps[column]!=0]
	return pdfeature





# Question: How many UCEs are there in a domain

# Question: Where in the domain are the UCEs

# Question: What other features characterize domains with UCEs










# map to features
def coverage_primary_to_secondary(btfeature,queryfeature):
	featurecov = queryfeature.coverage(btfeature)
	return featurecov

# get closest
def primary_to_tertiary_distances(btfeature,queryfeature):
	featuredis = btfeature.closest(queryfeature,d=True)
	return featuredis


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
	
	# 2) descriptive stats per feature file
	for file in allfiles:
		btfeature = get_bedtools_features(file)
		pdfeature = convert_bedtools_to_panda(btfeature)
		pdcoord = label_coordinate_columns(pdfeature)
		pdstats = descriptive_stats(pdcoord,'size')
	
	primary = get_bedtools_features(primaryfile)
	
	# 3) get overlaps - secondary
	for file in secondaryfiles:
		secondary = get_bedtools_features(file)
		intersect = intersect_bedfiles_c_true(secondary,primary)
		pdfeature = convert_bedtools_to_panda(intersect)
		pdcoord = label_coordinate_columns(pdfeature)
		pdcoord.columns.values[3]='intersect'
		count_number_with_zero_overlaps(pdcoord,'intersect')
		pdclean = remove_rows_with_no_overlaps(pdcoord,'intersect')
		intersectstats = descriptive_stats(pdclean,'intersect')
		print intersectstats
		
	
	# 4) get nearest - tertiary



if __name__ == "__main__":
	main()