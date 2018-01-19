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
	return pbt.BedTool(strFileName)

# convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	return pd.read_table(btfeature.fn,header=None)

# 4 cols coord labels
def label_coordinate_columns(pdfeature):
	pdfeature['size'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	return pdfeature

# get stats
def panda_describe_column(pdfeature,column):
	return pdfeature[column].describe()

# save panda
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	return afile.intersect(bfile,c=True)

# total number of elements without overlaps
def count_number_with_zero_overlaps(df,column):
	return len(df[(df[column]==0)])

# move elements without any overlaps
def remove_rows_with_no_overlaps(overlaps,column):
	return overlaps[overlaps[column]!=0]

# Query: Return some info about the element size stats from file
def descriptive_stats(file):
	btfeature = get_bedtools_features(file)
	pdfeature = convert_bedtools_to_panda(btfeature)
	pdcoord = label_coordinate_columns(pdfeature)
	pdstats = panda_describe_column(pdcoord,'size')
	print 'feature stats on element sizes for {0}'.format(file)
# 	print pdstats
	return pdstats

# Query: How many UCEs are there in a domain
def overlaping_features(primary,sfile):
	secondary = get_bedtools_features(sfile)
	intersect = intersect_bedfiles_c_true(secondary,primary)
	pdfeature = convert_bedtools_to_panda(intersect)
	pdcoord = label_coordinate_columns(pdfeature)
	pdcoord.columns.values[3]='intersect'
	nooverlaps = count_number_with_zero_overlaps(pdcoord,'intersect')
	print '{0} instances of no overlaps on {1}'.format(nooverlaps,sfile)
	pdclean = remove_rows_with_no_overlaps(pdcoord,'intersect')
	return secondary,pdclean

# intersect files and get original coords and overlap size
def intersect_bedfiles_wo_true(afile,bfile):
	return afile.intersect(bfile,wo=True)

# 7 cols coord labels
def label_expanded_coordinate_columns(pdfeature):
	pdfeature.columns.values[0]='achr'
	pdfeature.columns.values[1]='astart'
	pdfeature.columns.values[2]='aend'
	pdfeature.columns.values[3]='bchr'
	pdfeature.columns.values[4]='bstart'
	pdfeature.columns.values[5]='bend'
	pdfeature.columns.values[6]='overlapsize'
	return pdfeature

# make columns for distance between element and each boundary
def get_boundary_distances(pdfeature):
	pdfeature['enddistance'] = pdfeature['bend']-pdfeature['aend']
	pdfeature['startdistance'] = pdfeature['bstart']-pdfeature['astart']
	return pdfeature

# determine which boundary distance is further
def binary_boundary_distances(pdfeature):
	pdfeature['closerboundary'] = np.where((pdfeature['startdistance'] >= pdfeature['enddistance']),'start','end')
	return pdfeature

# Query: Where in the domain are the UCEs
def locateing_features(primary,file):
	intersect = intersect_bedfiles_wo_true(primary,file)
	pdfeature = convert_bedtools_to_panda(intersect)
	pdcoord = label_expanded_coordinate_columns(pdfeature)
	pdboundary = get_boundary_distances(pdcoord)
	pdbinary = binary_boundary_distances(pdboundary) # may want to bin domain instead, and return which bin the elements are in

# Query: What other features characterize domains with UCEs
def additional_features(secondary,sclean,tfile):
	tertiary = get_bedtools_features(tfile)
	intersect = intersect_bedfiles_c_true(secondary,tertiary)
	pdfeature = convert_bedtools_to_panda(intersect)
	pdcoord = label_coordinate_columns(pdfeature)
	pdcoord.columns.values[3]='intersect'
	print sclean
	print pdcoord.head()

def main():
	args = get_args()
	
	# read in files from args
	primaryfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tertiaryfiles = [line.strip() for line in args.tertiaryfeatures]
	
	genomefile = args.genomefile
	
	allfiles = secondaryfiles + tertiaryfiles
	allfiles.insert(0,primaryfile)
	
	print 'running script with files: {0}'.format(allfiles)
	
	# run descriptive stats per feature file
	for file in allfiles:
		pdstats = descriptive_stats(file)
	
	primary = get_bedtools_features(primaryfile)
	
	# process feature files
	for sfile in secondaryfiles:
		secondary,sclean = overlaping_features(primary,sfile)
		print 'feature overlap stats on {0}'.format(sfile)
		intersectstats = panda_describe_column(sclean,'intersect')
# 		print intersectstats
		locateing_features(primary,sfile)
		for tfile in tertiaryfiles:
			additional_features(secondary,sclean,tfile)








if __name__ == "__main__":
	main()