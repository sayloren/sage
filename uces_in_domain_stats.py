"""
Script to print file about number of uces and genes within a domain stats

Wren Saylor
Feb 2018

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

In:
primary - uce file
secondary - file with list of domain filenames
tertiary - genes

Out:
stats file with the intersects; uce x domain, gene x domain, domain size
stats file for all the domains; uce size, all domain size, gene size

"""
import argparse
import pandas as pd
import pybedtools as pbt
from itertools import cycle
import numpy as np

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",required=True,type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeature",type=str,help="the tertiary elements file")# Genes
	return parser.parse_args()

# get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

# convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	return pd.read_table(btfeature.fn,header=None)

# print the number of elements in the feature file
def print_number_of_elements(lengthfeatures,file):
	print 'there are {0} elements in {1}'.format(lengthfeatures,file)

# intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	return afile.intersect(bfile,c=True)

def run_print_number_file_features(file):
	btfeatures = get_bedtools_features(file)
	pdfeatures = convert_bedtools_to_panda(btfeatures)
	pdlabel = label_coordinate_columns(pdfeatures)
	print_number_of_elements(len(pdfeatures),file)
	return pdlabel

# coordinate labels and size
def label_coordinate_columns(pdfeature):
	pdfeature['size'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	return pdfeature

# create panda for overlap count datasets
def count_overlap_df(secondary,file,label):
	pdintersect = intersect_bedfiles_c_true(secondary,file)
	pdfeatures = convert_bedtools_to_panda(pdintersect)
	pdcoordinates = label_coordinate_columns(pdfeatures)
	pdcoordinates.columns.values[3]='intersect_{0}'.format(label)
	return pdcoordinates

# total number of elements without overlaps
def count_number_with_zero_overlaps(df,column):
	return len(df[(df[column]==0)])

# move elements without any overlaps
def remove_rows_with_no_overlaps(overlaps,column):
	return overlaps[overlaps[column]!=0]

# get stats for list of columns
def panda_describe_multiple_column(pdfeature):
	intersectcols = [col for col in pdfeature.columns if 'intersect' in col]
	sizecol = [col for col in pdfeature.columns if 'size' in col]
	statcols=intersectcols+sizecol
	return pdfeature[statcols].describe()

# save panda to file with mode 'a' for appending
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True,mode='a')

# get stats for single column
def panda_describe_single_column(pdfeatures,name):
	return pdfeatures[name].describe()

def main():
	args = get_args()
	
	# read in files from args
	pfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tfile = args.tertiaryfeature
	
	# print the number of features in single feature files
	labelprimary = run_print_number_file_features(pfile)
	labeltertiary = run_print_number_file_features(tfile)

	# initiate collection
	lumpsecondary = []
	
	# process feature files
	for sfile in secondaryfiles:
	
		# print number of features in secondary file
		run_print_number_file_features(sfile)
		
		# get secondary features
		secondary = get_bedtools_features(sfile)
		
		# make the pandas data sets for the count overlaps
		pdprimary = count_overlap_df(secondary,pfile,'{0}'.format(pfile))
		pdtertiary = count_overlap_df(secondary,tfile,'{0}'.format(tfile))
		
		# print the number of domains that do not have any uces
		nooverlaps = count_number_with_zero_overlaps(pdprimary,'intersect_{0}'.format(pfile))
		print '{0} instances of no overlaps of primary element on {1}'.format(nooverlaps,sfile)
		
		# concat the three data sets together
		concattotal = pdprimary.merge(pdtertiary,how='inner',on=['chr','start','end','size'])
		
		# remove the regions with no uces
		cleantotal = remove_rows_with_no_overlaps(concattotal,'intersect_{0}'.format(pfile))
		
		# add the data set to the lump sum to get the total stats at the end of the script
		lumpsecondary.append(cleantotal)
		
		# describe any file with a partial match of 'size' or 'intersect'
		cleandescribe = panda_describe_multiple_column(cleantotal)
		
		# add the domain file name to the columns
		cleandescribe.columns = [str(col) + '_{0}'.format(sfile) for col in cleandescribe.columns]
		
		# save panda to file
		save_panda(cleandescribe,'stats_{0}_individual_domains.txt'.format(pfile))
		
	# get the stats for the individual files
	primarystats = panda_describe_single_column(labelprimary,'size')
	tertiarystats = panda_describe_single_column(labeltertiary,'size')
	
	# concat and get stats for all the secondary files
	secondaryconcat = pd.concat(lumpsecondary)
	secondarystats = panda_describe_single_column(secondaryconcat,'size')
	
	# get the total count stats for all the secondary files
	primarycountstats = panda_describe_single_column(secondaryconcat,'intersect_{0}'.format(pfile))
	tertiarycountstats = panda_describe_single_column(secondaryconcat,'intersect_{0}'.format(tfile))
	
	# concat all the stats into panda
	allstats = pd.concat([primarystats,secondarystats,tertiarystats,primarycountstats,tertiarycountstats],axis=1)
	
	# correct column names with file names for clarity
	allstats.columns = ['size_{0}'.format(pfile),'size_all_domains','size_{0}'.format(tfile),'counts_in_UCE_domains_{0}'.format(pfile),'counts_in_UCE_domains_{0}'.format(tfile)]
	
	# save the panda to file
	save_panda(allstats,'stats_{0}_all_domains.txt'.format(pfile))

if __name__ == "__main__":
	main()