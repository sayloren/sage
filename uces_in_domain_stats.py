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
from scipy import stats

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
	print_number_of_elements(len(pdfeatures),file)
	return label_coordinate_columns(pdfeatures)

# coordinate labels and size
def label_coordinate_columns(pdfeature):
	pdfeature['Size(Kb)'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	pdfeature['Size(Kb)'] /= 1000.# convert to Kb
	return pdfeature

# create panda for overlap count datasets
def count_overlap_df(secondary,file,label):
	pdintersect = intersect_bedfiles_c_true(secondary,file)
	pdfeatures = convert_bedtools_to_panda(pdintersect)
	pdcoordinates = label_coordinate_columns(pdfeatures)
	pdcoordinates.columns.values[3]='intersect_{0}'.format(label)
	return pdcoordinates

# split the with and without element groups
def split_rows_with_no_overlaps(overlaps,column):
	withuces = overlaps[overlaps[column]!=0]
	withoutuces = overlaps[overlaps[column]==0]
	return withuces,withoutuces

# get stats for list of columns
def panda_describe_multiple_column(pdfeature):
	intersectcols = [col for col in pdfeature.columns if 'intersect' in col]
	sizecol = [col for col in pdfeature.columns if 'Size(Kb)' in col]
	densitycol = [col for col in pdfeature.columns if 'Gene_Density(Kb)' in col]
	statcols=intersectcols+sizecol+densitycol
	return pdfeature[statcols].describe()

# select a list of columns where the column name contains a key word
def subset_column_by_keyword(pdfeature,key):
	return pdfeature.filter(regex=key)

# save panda to file with mode 'a' for appending
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True,mode='a')#

# get stats for single column
def panda_describe_single_column(pdfeatures,name):
	return pdfeatures[name].describe()

# split into the three stats files and save
def split_key_words_to_files(pdfeatures,pdfeaturesout,pval,pfile):
	# domain size
	sizedescribe = subset_column_by_keyword(pdfeatures,'Size')
	sizedescribeout = subset_column_by_keyword(pdfeaturesout,'Size')
	sizedescribepval = subset_column_by_keyword(pval,'Size')
	save_panda(sizedescribe,'stats_domain_size_{0}.txt'.format(pfile))
	save_panda(sizedescribeout,'stats_domain_size_{0}.txt'.format(pfile))
	save_panda(sizedescribepval,'stats_domain_size_{0}.txt'.format(pfile))
	# gene density
	densitydescribe = subset_column_by_keyword(pdfeatures,'Gene_Density')
	densitydescribeout = subset_column_by_keyword(pdfeaturesout,'Gene_Density')
	densitydescribepval = subset_column_by_keyword(pval,'Gene_Density')
	save_panda(densitydescribe,'stats_gene_density_{0}.txt'.format(pfile))
	save_panda(densitydescribeout,'stats_gene_density_{0}.txt'.format(pfile))
	save_panda(densitydescribepval,'stats_gene_density_{0}.txt'.format(pfile))
	# number UCEs
	intersectdescribe = subset_column_by_keyword(pdfeatures,'intersect')
	intersectdescribeout = subset_column_by_keyword(pdfeaturesout,'intersect')
	intersectdescribepval = subset_column_by_keyword(pval,'intersect')
	save_panda(intersectdescribe,'stats_number_uces_{0}.txt'.format(pfile))
	save_panda(intersectdescribeout,'stats_number_uces_{0}.txt'.format(pfile))
	save_panda(intersectdescribepval,'stats_number_uces_{0}.txt'.format(pfile))

def main():
	args = get_args()
	
	pfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tfile = args.tertiaryfeature
	
	labelprimary = run_print_number_file_features(pfile)
	labeltertiary = run_print_number_file_features(tfile)
	
	lumpsec,lumpsecstats = [],[]
	lumpsecout,lumpsecstatsout = [],[]
	lumppval = []
	
	for sfile in secondaryfiles:
		
		run_print_number_file_features(sfile)
		secondary = get_bedtools_features(sfile)
		pdprimary = count_overlap_df(secondary,pfile,'UCEs')
		pdtertiary = count_overlap_df(secondary,tfile,'Genes')
		
		concattotal = pdprimary.join(pdtertiary,rsuffix='_extra')
		concattotal['Gene_Density(Kb)'] = concattotal['intersect_Genes']/concattotal['Size(Kb)']
		concattotal.drop(columns=['chr_extra','start_extra','end_extra','Size(Kb)_extra','intersect_Genes'],axis=1,inplace=True)
		
		withuces,withoutuces = split_rows_with_no_overlaps(concattotal,'intersect_UCEs')
		
		statcoefintersect,statpvalintersect = stats.mannwhitneyu(withuces['intersect_UCEs'],withoutuces['intersect_UCEs'])
		statcoefsize,statpvalsize = stats.mannwhitneyu(withuces['Size(Kb)'],withoutuces['Size(Kb)'])
		statcoefdensity,statpvaldensity = stats.mannwhitneyu(withuces['Gene_Density(Kb)'],withoutuces['Gene_Density(Kb)'])
		formatpvalintersect = '{:.02e}'.format(statpvalintersect)
		formatpvalsize = '{:.02e}'.format(statpvalsize)
		formatpvaldensity = '{:.02e}'.format(statpvaldensity)
		
		pdmannwhitney = pd.DataFrame({'p-value':[formatpvalintersect,formatpvalsize,formatpvaldensity],'coeff':[statcoefintersect,statcoefsize,statcoefdensity],'dataset':['intersect_UCEs_{0}'.format(sfile),'Size(Kb)_{0}'.format(sfile),'Gene_Density(Kb)_{0}'.format(sfile)]})
		lumppval.append(pdmannwhitney)
		
		lumpsec.append(withuces)
		lumpsecout.append(withoutuces)
		
		withdescribe = panda_describe_multiple_column(withuces)
		withoutdescribe = panda_describe_multiple_column(withoutuces)
		
		withdescribe.columns = [str(col) + '_{0}'.format(sfile) for col in withdescribe.columns]
		withoutdescribe.columns = [str(col) + '_without_{0}'.format(sfile) for col in withoutdescribe.columns]
		
		lumpsecstats.append(withdescribe)
		lumpsecstatsout.append(withoutdescribe)
	
	seccat = pd.concat(lumpsec)
	seccatout = pd.concat(lumpsecout)
	
	statcoefintersect,statpvalintersect = stats.mannwhitneyu(seccat['intersect_UCEs'],seccatout['intersect_UCEs'])
	statcoefsize,statpvalsize = stats.mannwhitneyu(seccat['Size(Kb)'],seccatout['Size(Kb)'])
	statcoefdensity,statpvaldensity = stats.mannwhitneyu(seccat['Gene_Density(Kb)'],seccatout['Gene_Density(Kb)'])
	formatpvalintersect = '{:.02e}'.format(statpvalintersect)
	formatpvalsize = '{:.02e}'.format(statpvalsize)
	formatpvaldensity = '{:.02e}'.format(statpvaldensity)
	
	pdmannwhitney = pd.DataFrame({'p-value':[formatpvalintersect,formatpvalsize,formatpvaldensity],'coeff':[statcoefintersect,statcoefsize,statcoefdensity],'dataset':['intersect_UCEs_All_Domains','Size(Kb)_All_Domains','Gene_Density(Kb)_All_Domains']})
	lumppval.append(pdmannwhitney)
	
	secdescribe = panda_describe_multiple_column(seccat)
	secdescribeout = panda_describe_multiple_column(seccatout)
	
	secdescribe.columns = [str(col) + '_All_Domains' for col in secdescribe.columns]
	secdescribeout.columns = [str(col) + '_without_All_Domains' for col in secdescribeout.columns]
	
	lumpsecstats.append(secdescribe)
	lumpsecstatsout.append(secdescribeout)
	
# 	labelprimary.rename(columns={'Size(Kb)':'Size(Kb)_{0}'.format(pfile)},inplace=True)
# 	primarystats = panda_describe_single_column(labelprimary,'Size(Kb)_{0}'.format(pfile))
# 	labeltertiary.rename(columns={'Size(Kb)':'Size(Kb)_{0}'.format(tfile)},inplace=True)
# 	tertiarystats = panda_describe_single_column(labeltertiary,'Size(Kb)_{0}'.format(tfile))
# 	lumpsecstats.append(primarystats)
# 	lumpsecstats.append(tertiarystats)
	secstatscat = pd.concat(lumpsecstats,axis=1)
	secstatscatout = pd.concat(lumpsecstatsout,axis=1)
	pvalcat = pd.concat(lumppval)
	pvalcat.set_index(keys='dataset',inplace=True,drop=True)
	tpvalcat = pvalcat.T
	
	split_key_words_to_files(secstatscat,secstatscatout,tpvalcat,pfile)

if __name__ == "__main__":
	main()