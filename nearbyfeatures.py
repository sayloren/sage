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

Outputs:
File with stats about (by appearance in column labels):
	-intersect between primary overlaping with secondary
	-intersect between tertiary overlapping secondary
	-size of the secondary regions
	-size of the primary region
	-size of the tertiary region
Bar plot for each secondary file grouped with primary and tertiary files by which bin in the seconary they fall

"""
import argparse
import pandas as pd
import pybedtools as pbt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import cycle
import matplotlib

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",required=True,type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the tertiary features to query")# Genes
	parser.add_argument("-q","--quinaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the quinary features to query")# Mouse UCEs
	parser.add_argument("-g","--genomefile",type=str,help="genome file",default='hg19.genome')
	parser.add_argument("-b","--binnumber",type=int,default='10',help='number of bins to chunk the secondary files into, must be even number')
	return parser.parse_args()

# 1a,3a) get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

# 2a) convert panda to bedtool
def convert_panda_to_bed_format(panda):
	arArFeatures = panda.values.tolist()
	btoutFeatures = pbt.BedTool(arArFeatures)
	return btoutFeatures

# 2b) bin secondary regions
def make_window_with_secondary_files(sfile,bins):
	windows = pbt.BedTool().window_maker(b=sfile,n=bins,i="src")
	return windows

# 1b,2c,3b,4a) intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	return afile.intersect(bfile,c=True)

# 4b) intersect files and get original coords and overlap size
def intersect_bedfiles_wo_true(afile,bfile):
	return afile.intersect(bfile,wo=True)

# 1c,2d,3c,4c) convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	return pd.read_table(btfeature.fn,header=None)

# 1d,2e,3d) 4 cols coord labels
def label_coordinate_columns(pdfeature):
	pdfeature['size'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	return pdfeature

# 2f) merge a and b files on 'id' column
def intersect_pandas_with_id(afile,bfile):
	return pd.merge(afile,bfile,on='id')

# 1e) save pandas in bedtools format
def save_panda_bed_format(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=False,header=False)

# 2f,4d) 7 cols coord labels
def label_expanded_coordinate_columns(pdfeature):
	pdfeature.columns = ['achr','astart','aend','countfeature','bchr','bstart','bend','overlapsize']
	return pdfeature

# 4e) groupby larger region
def group_df_by_secondary_regions(pdfeature):
	countcolumns = [col for col in pdfeature.columns if 'count' in col]
	outgroup = []
	for col in countcolumns:
		group = pd.DataFrame({'bincounts':pdfeature.groupby(['bchr','bstart','bend'])[col].apply(list)}).reset_index()
		group.columns = ['chr','start','end','bincounts_{0}'.format(col)]
		outgroup.append(group)
	grouped = reduce(lambda x, y: pd.merge(x,y,on=['chr','start','end']),outgroup)
	return grouped

# 5a) get stats for single column
def panda_describe_single_column(btfeature,name):
	pdfeature = convert_bedtools_to_panda(btfeature)
	pdlabel = label_coordinate_columns(pdfeature)
	pdselect = pdlabel[['chr','start','end','size']]
	pdselect.rename(columns={'size':'{0}'.format(name)},inplace=True)
	statcol = [col for col in pdselect.columns if 'size' in col]
	describe = pdselect[statcol].describe()
	return describe

# 5b) total number of elements without overlaps
def count_number_with_zero_overlaps(df,column):
	return len(df[(df[column]==0)])

# 5c) move elements without any overlaps
def remove_rows_with_no_overlaps(overlaps,column):
	return overlaps[overlaps[column]!=0]

# 5d) get stats for list of columns
def panda_describe_multiple_column(pdfeature):
	intersectcols = [col for col in pdfeature.columns if 'intersect' in col]
	sizecol = [col for col in pdfeature.columns if 'size' in col]
	statcols=intersectcols+sizecol
	return pdfeature[statcols].describe()

# 5e) save panda
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True,mode='a')

# 1) Query: How many UCEs are there in a domain 
def overlaping_features(primary,sfile):
	secondary = get_bedtools_features(sfile)
	intersect = intersect_bedfiles_c_true(secondary,primary)
	pdfeature = convert_bedtools_to_panda(intersect)
	pdcoord = label_coordinate_columns(pdfeature)
	pdcoord.columns.values[3]='intersect_primary'
	pdcoord.insert(len(pdcoord.columns),'id',range(0,0+len(pdcoord)))
	# save scoord as temporary file to read in later with id
	return secondary,pdcoord

# 2) Query: Where in the domain are the UCEs - going to want column 'count_primary_overlaps_windows'
def locateing_features(primary,sfile,bins):
	sfeature = sfile[['chr','start','end','id']] # get just the coordinates and the id
	sbedtool = convert_panda_to_bed_format(sfeature) # convert to bedtool
	windows = make_window_with_secondary_files(sbedtool,bins) # make the windows with the secondary file, perserved the id
	intersectprimary = intersect_bedfiles_c_true(windows,primary) # get the number of primary elements at each window
	
	primarypd = convert_bedtools_to_panda(intersectprimary)
	primarylabel = label_coordinate_columns(primarypd)
	primarylabel.columns.values[3]='id'
	primarylabel.columns.values[4]='count_primary_overlaps_windows'
	primarylabel.rename(columns={'chr':'window_chr','start':'window_start','end':'window_end','size':'window_size'},inplace=True)
	sfile.rename(columns={'chr':'secondary_chr','start':'secondary_start','end':'secondary_end','size':'secondary_size','intersect_primary':'count_primary_overlaps_secondary'},inplace=True)

	# this is where rao is having trouble: the domains are redundant - have to have an id to merge on
	# old method commented out
# 	intersectsecondary = intersect_bedfiles_wo_true(intersectprimary,sfile)
# 	pdfeature = convert_bedtools_to_panda(intersectsecondary)
# 	pdcoord = label_expanded_coordinate_columns(pdfeature)
# 	pdcoord.rename(columns={'countfeature':'count_primary'},inplace=True)
	
	intersectsecondary = intersect_pandas_with_id(primarylabel,sfile)
	
	return intersectsecondary,windows

# 3) Query: What other features characterize domains with UCEs; in a domain
def additional_features_overlap(secondary,scoord,tfile):
	tertiary = get_bedtools_features(tfile)
	intersect = intersect_bedfiles_c_true(secondary,tertiary)
	pdfeature = convert_bedtools_to_panda(intersect)
	tcoord = label_coordinate_columns(pdfeature)
	tcoord.columns.values[3]='intersect_{0}'.format(tfile)
	concatintersect = pd.merge(scoord,tcoord,how='inner',on=['chr','start','end'])
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
	return pdgroup

# 6a) remove the rows where there are no primary element overlaps with the secondary regions
def drop_primary_zero_list(pdfeatures,column):
	thresh = pdfeatures[~pdfeatures[column].apply(lambda row: all(item ==0 for item in row))]
	return thresh

# 6b) keep only those rows where there are no primary element overlaps with the secondary regions
def keep_primary_zero_list(pdfeatures,column):
	thresh = pdfeatures[pdfeatures[column].apply(lambda row: all(item ==0 for item in row))]
	return thresh

# 6c) format the binned data frame for easy graphing
def format_binned_data_for_graphing(pdfeatures,bins):
	selectcols = [col for col in pdfeatures.columns if 'bincounts' in col]
	list = []
	for group in selectcols:
		subset = pdfeatures[[group]]
		split = pd.DataFrame(subset[group].values.tolist())
		sum = split.sum(axis=0)
		list.append(sum)
	concat = pd.concat(list,axis=1)
	concat.columns = [selectcols]
	bincolumns = range(bins)
	format = pd.melt(concat)
	format.columns = ['filename','sumbin']
	format['bin'] = bincolumns * (format.shape[0]/len(bincolumns))
	return format

# 6d) fold data set
def fold_formated_binned_data(pdfeatures,bins):
	halfbin = bins/2
	pdfeatures['invertsumbin'] = pdfeatures['sumbin'].iloc[::-1].reset_index(drop=True)
	pdfeatures['sumsums'] = pdfeatures['invertsumbin'] + pdfeatures['sumbin']
	headfeatures = pdfeatures.head(n=halfbin)
	dropfeatures = headfeatures[['filename','sumsums','bin']]
	dropfeatures.columns = ['filename','sumbin','bin']
	return dropfeatures

# 6ei) graph box plots for secondary and tertiary data - will have to put in the y axis label as arg
def graph_boxplot_region_size(pdfeatures,primaryfile,sfile,tfile,ylabeltext):
	sns.set_style('ticks')
	pp = PdfPages('tertiarycount_{0}_{1}_{2}.pdf'.format(primaryfile,sfile,tfile))
	plt.figure(figsize=(14,7))
	
	sns.set_palette("Blues")
	
	gs = gridspec.GridSpec(1,1)
	gs.update(hspace=.8)
	
	ax0 = plt.subplot(gs[0,0])
	sns.boxplot(data=pdfeatures,x='primary',y=tfile,showfliers=False)
	ax0.set_ylabel(ylabeltext)
	ax0.set_xlabel('')
# 	ax0.tick_params(axis='both',which='major',labelsize=16)
	for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] + ax0.get_xticklabels() + ax0.get_yticklabels()):
		item.set_fontsize(22)
	
	sns.despine()
	pp.savefig()
	pp.close()

# 6eii) line graph binned data with no tertiary features
def graph_binned_regions_no_tertiary(pdfeatures,primaryfile,sfile):
	sns.set_style('ticks')
	pp = PdfPages('bincounts_{0}_{1}.pdf'.format(primaryfile,sfile))
	plt.figure(figsize=(14,7))
	
	unique = len(pdfeatures['filename'].unique())
	sns.set_palette("Blues",n_colors=unique)
	
	gs = gridspec.GridSpec(1,1)
	gs.update(hspace=.8)
	
	ax0 = plt.subplot(gs[0,0])
	sns.pointplot(data=pdfeatures,x='bin',y='sumbin',color='#9ecae1',scale=3)
	ax0.set_ylabel('Frequency')
	ax0.set_xlabel('Bin Distance from Edge')
# 	ax0.tick_params(axis='both',which='major',labelsize=16)
	for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] + ax0.get_xticklabels() + ax0.get_yticklabels()):
		item.set_fontsize(22)
	
	sns.despine()
	pp.savefig()
	pp.close()

def main():
	args = get_args()
	
	# read in files from args
	primaryfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tertiaryfiles = [line.strip() for line in args.tertiaryfeatures]
	quinaryfiles = [line.strip() for line in args.quinaryfeatures] # have to integrate to graph subset of uces on boxplot
	bins = args.binnumber
	genomefile = args.genomefile
	
	# make list of all files to iterate through
	allfiles = secondaryfiles + tertiaryfiles
	allfiles.insert(0,primaryfile)
	
	primary = get_bedtools_features(primaryfile)
	
	lumpsecondaryfiles = []
	# process feature files
	for sfile in secondaryfiles:
		# 1) Query: How many UCEs are there in a domain
		secondary,scoord = overlaping_features(primary,sfile)
		
		# get the number of secondary features with no primary overlaps
		nooverlaps = count_number_with_zero_overlaps(scoord,'intersect_primary')
		print '{0} instances of no overlaps of primary element on {1}'.format(nooverlaps,sfile)
		
		# remove all secondary regions with no primary overlaps
		cleanintersect = remove_rows_with_no_overlaps(scoord,'intersect_primary') # may want to wait, and make graphs for those with uces vs those without
		
		# 2a) Query: Where in the domain are the UCEs
		binfeatures,windows = locateing_features(primary,scoord,bins)
		
# 		# group bin counts by secondary features they came from
		groupfeatures = group_df_by_secondary_regions(binfeatures)
# 		
		# remove secondary features where there are no primary elements
# 		dropzero = drop_primary_zero_list(groupfeatures,'bincounts_count_primary')
# 		
# 		# shape the data frame
# 		formatbin = format_binned_data_for_graphing(dropzero,bins)
# 		
# 		# fold the data frame by combining edges
# 		foldbin = fold_formated_binned_data(formatbin,bins)
# 		
# 		# 6) generate graph for bin results
# 		graph_binned_regions_no_tertiary(foldbin,primaryfile,sfile)
# 		
# 		# format data for boxplot graphs of secondary size
# 		dropzero = (scoord.loc[scoord['intersect_primary'] > 0])
# 		sumdrop = dropzero['size'].reset_index(drop=False)
# 		sumdrop['primary'] = 'Regions With Primary'
# 		keepzero = (scoord.loc[scoord['intersect_primary'] < 1])
# 		sumkeep = keepzero['size'].reset_index(drop=False)
# 		sumkeep['primary'] = 'Regions Without Primary'
# 		tertiarysum = pd.concat([sumdrop,sumkeep])
# 		tertiarysum.drop(columns=['index'],inplace=True)
# 		
# 		# graph size for secondary regions
# 		graph_boxplot_region_size(tertiarysum,primaryfile,sfile,'size','Size (kb)')
# 		
# 		# add the primary x secondary intersections data to the list for lump secondary stats
# 		lumpsecondaryfiles.append(cleanintersect)
# 		
		# for contrast with tertiary features
# 		for tfile in tertiaryfiles:
# 			# 3) Query: What other features characterize domains with UCEs; in a domain
# 			tertiary,intersect = additional_features_overlap(secondary,scoord,tfile)
# 			
# 			# 4) Query: What other features characterize domains with UCEs; where in domain
# 			groupfeatures = additional_features_locate(secondary,tertiary,tfile,binfeatures,windows)
# 			
# 			# 5) generate stats results
# 			primarystats = panda_describe_single_column(primary,'size_primary')
# 			tertiarystats = panda_describe_single_column(tertiary,'size_{0}'.format(tfile))
# 			
# 			# get the number of features with no primary overlaps
# 			nooverlaps = count_number_with_zero_overlaps(intersect,'intersect_primary')
# 			print '{0} instances of no overlaps of primary element on {1}'.format(nooverlaps,sfile)
# 			
# 			# remove all secondary regions with no primary overlaps
# 			cleanintersect = remove_rows_with_no_overlaps(intersect,'intersect_primary') # may want to wait, and make graphs for those with uces vs those without
# 			cleanintersect.rename(columns={'size_x':'size_{0}'.format(sfile)},inplace=True)
# 			cleanintersect.drop(columns=['size_y'],inplace=True)
# 			
# 			# make stats file
# 			intersectstats = panda_describe_multiple_column(cleanintersect)
# 			alltertiarystats = pd.concat([intersectstats,primarystats,tertiarystats],axis=1)
# 			
# 			# save stats to file
# 			save_panda(alltertiarystats,'stats_{0}_intersection.txt'.format(primaryfile))
# 			
# 			# format data for boxplot graphs of tertiary feature counts
# 			dropzero = (intersect.loc[intersect['intersect_primary'] > 0])
# 			sumdrop = dropzero['intersect_{0}'.format(tfile)].reset_index(drop=False)
# 			sumdrop['primary'] = 'Regions With Primary'
# 			keepzero = (intersect.loc[intersect['intersect_primary'] < 1])
# 			sumkeep = keepzero['intersect_{0}'.format(tfile)].reset_index(drop=False)
# 			sumkeep['primary'] = 'Regions Without Primary'
# 			tertiarysum = pd.concat([sumdrop,sumkeep])
# 			tertiarysum.drop(columns=['index'],inplace=True)
# 			
# 			# graph counts for tertiary features
# 			graph_boxplot_region_size(tertiarysum,primaryfile,sfile,'intersect_{0}'.format(tfile),'Frequency')
# 			
# 	# make the stats data frame for all the secondary features
# 	concatsecondary = pd.concat(lumpsecondaryfiles)
# 	concatsecondary.rename(columns={'size':'size_all_secondary'},inplace=True)
# 	catsecondarystats = panda_describe_multiple_column(concatsecondary)
# 	
# 	# save stats to file
# 	save_panda(catsecondarystats,'stats_{0}_intersection.txt'.format(primaryfile))

if __name__ == "__main__":
	main()