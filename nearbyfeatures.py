"""

Script to look at how close features are between two bed files, or a UCE file and a list of bedfiles

Wren Saylor
Jan 2018

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

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the tertiary features to query")# Genes
	parser.add_argument("-g","--genomefile",type=str,help="genome file",default='hg19.genome')
	parser.add_argument("-b","--binnumber",type=int,default='10',help='number of bins to chunk the secondary files into')
	return parser.parse_args()

# 1a,3a) get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

# 2a) bin secondary regions
def make_window_with_secondary_files(sfile,bins):
	windows = pbt.BedTool().window_maker(b=sfile,n=bins)
	return windows

# 1b,2b,3b,4a) intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	return afile.intersect(bfile,c=True)

# 2c,4b) intersect files and get original coords and overlap size
def intersect_bedfiles_wo_true(afile,bfile):
	return afile.intersect(bfile,wo=True)

# 1c,2d,3c,4c) convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	return pd.read_table(btfeature.fn,header=None)

# 1d,3d) 4 cols coord labels
def label_coordinate_columns(pdfeature):
	pdfeature['size'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	return pdfeature

# 2e,4d) 7 cols coord labels
def label_expanded_coordinate_columns(pdfeature):
	pdfeature.columns = ['achr','astart','aend','countfeature','bchr','bstart','bend','overlapsize']
	return pdfeature

# groupby larger region
def group_df_by_secondary_regions(pdfeature):
# 	grouped = pdfeature.groupby(['bchr','bstart','bend'])['count_primary'].apply(list)
	countcolumns = [col for col in pdfeature.columns if 'count' in col]
	outgroup = []
	for col in countcolumns:
		group = pdfeature.groupby(['bchr','bstart','bend'])[col].apply(list)
		outgroup.append(group)
	grouped = pd.concat(outgroup)
	print grouped
	return grouped




# get stats
def panda_describe_column(pdfeature,column):
	return pdfeature[column].describe()

# total number of elements without overlaps
def count_number_with_zero_overlaps(df,column):
	return len(df[(df[column]==0)])

# move elements without any overlaps
def remove_rows_with_no_overlaps(overlaps,column):
	return overlaps[overlaps[column]!=0]

# save panda
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True)

# 1) Query: How many UCEs are there in a domain
def overlaping_features(primary,sfile):
	secondary = get_bedtools_features(sfile)
	intersect = intersect_bedfiles_c_true(secondary,primary)
	pdfeature = convert_bedtools_to_panda(intersect)
	pdcoord = label_coordinate_columns(pdfeature)
	pdcoord.columns.values[3]='intersect_primary'
	return secondary,pdcoord

# 2) Query: Where in the domain are the UCEs
def locateing_features(primary,sfile,bins):
	windows = make_window_with_secondary_files(sfile,bins)
	intersectprimary = intersect_bedfiles_c_true(windows,primary)
	intersectsecondary = intersect_bedfiles_wo_true(intersectprimary,sfile)
	pdfeature = convert_bedtools_to_panda(intersectsecondary)
	pdcoord = label_expanded_coordinate_columns(pdfeature)
	pdcoord.rename(columns={'countfeature':'count_primary'},inplace=True)
# 	pdgroup = group_df_by_secondary_regions(pdcoord)
	return pdcoord,windows

# 3) Query: What other features characterize domains with UCEs; in a domain
def additional_features_overlap(secondary,scoord,tfile):
	tertiary = get_bedtools_features(tfile)
	intersect = intersect_bedfiles_c_true(secondary,tertiary)
	pdfeature = convert_bedtools_to_panda(intersect)
	tcoord = label_coordinate_columns(pdfeature)
	tcoord.columns.values[3]='intersect_{0}'.format(tfile)
	concatintersect = pd.merge(scoord,tcoord,how='inner',on=['chr','start','end','size'])
	return tertiary,concatintersect

# 4) Query: What other features characterize domains with UCEs; where in domain
def additional_features_locate(secondary,tertiary,tfile,scoord,windows):
	intersecttertiary = intersect_bedfiles_c_true(windows,tertiary)
	intersectsecondary = intersect_bedfiles_wo_true(intersecttertiary,secondary)
	pdfeature = convert_bedtools_to_panda(intersectsecondary)
	tcoord = label_expanded_coordinate_columns(pdfeature)
	tcoord.rename(columns={'countfeature':'count_{0}'.format(tfile)},inplace=True)
	concatintersect = pd.merge(scoord,tcoord,how='inner',on=['achr','astart','aend','bchr','bstart','bend'])
	pdgroup = group_df_by_secondary_regions(concatintersect)

def main():
	args = get_args()
	
	# read in files from args
	primaryfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tertiaryfiles = [line.strip() for line in args.tertiaryfeatures]
	bins = args.binnumber
	genomefile = args.genomefile
	
	# make list of all files to iterate through
	allfiles = secondaryfiles + tertiaryfiles
	allfiles.insert(0,primaryfile)
	
	primary = get_bedtools_features(primaryfile)
	
	# process feature files
	for sfile in secondaryfiles:
		secondary,scoord = overlaping_features(primary,sfile)
		binfeatures,windows = locateing_features(primary,sfile,bins)
		for tfile in tertiaryfiles:
			tertiary,intersect = additional_features_overlap(secondary,scoord,tfile)
			additional_features_locate(secondary,tertiary,tfile,binfeatures,windows)
#			return some stats from file; size, intersections by type
# 			intersectstats = panda_describe_column(coords,['size','intersect_primary','intersect...types'])
# 			nooverlaps = count_number_with_zero_overlaps(intersect,'intersect_primary')
# 			print '{0} instances of no overlaps on {1}'.format(nooverlaps,sfile)
# 			pdclean = remove_rows_with_no_overlaps(intersect,'intersect_primary') # may want to wait, and make graphs for those with uces vs those without
			

if __name__ == "__main__":
	main()