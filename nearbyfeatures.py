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
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import cycle

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
	return secondary,pdcoord

# 2) Query: Where in the domain are the UCEs
def locateing_features(primary,sfile,bins):
	windows = make_window_with_secondary_files(sfile,bins)
	intersectprimary = intersect_bedfiles_c_true(windows,primary)
	intersectsecondary = intersect_bedfiles_wo_true(intersectprimary,sfile)
	pdfeature = convert_bedtools_to_panda(intersectsecondary)
	pdcoord = label_expanded_coordinate_columns(pdfeature)
	pdcoord.rename(columns={'countfeature':'count_primary'},inplace=True)
	return pdcoord,windows

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
def drop_primary_zero_list(pdfeatures):
	thresh = pdfeatures[~pdfeatures['bincounts_count_primary'].apply(lambda row: all(item ==0 for item in row))]
	return thresh

# 6b) format the binned data frame for easy graphing
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

# 6d) graph binned data
def graph_binned_regions(pdfeatures,primaryfile,sfile):
	sns.set_style('ticks')
	pp = PdfPages('bincounts_{0}_{1}.pdf'.format(primaryfile,sfile))
	plt.figure(figsize=(14,7))
	
	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	gs.update(hspace=.8)
	
	ax0 = plt.subplot(gs[0,:])
	
	sns.barplot(data=pdfeatures,x='bin',y='sumbin',hue='filename')
	
	ax0.set_ylabel('Counts for {0}'.format(sfile),size=16)
	ax0.tick_params(axis='both',which='major',labelsize=16)

	sns.despine()
	pp.savefig()
	pp.close()

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
		# 1) Query: How many UCEs are there in a domain
		secondary,scoord = overlaping_features(primary,sfile)
		
		# 2) Query: Where in the domain are the UCEs
		binfeatures,windows = locateing_features(primary,sfile,bins)
		
		for tfile in tertiaryfiles:
			# 3) Query: What other features characterize domains with UCEs; in a domain
			tertiary,intersect = additional_features_overlap(secondary,scoord,tfile)
			
			# 4) Query: What other features characterize domains with UCEs; where in domain
			groupfeatures = additional_features_locate(secondary,tertiary,tfile,binfeatures,windows)
			
			# 5) generate stats results
			primarystats = panda_describe_single_column(primary,'size_primary')
			tertiarystats = panda_describe_single_column(tertiary,'size_{0}'.format(tfile))
			
			nooverlaps = count_number_with_zero_overlaps(intersect,'intersect_primary')
			print '{0} instances of no overlaps of primary element on {1}'.format(nooverlaps,sfile)
			
			cleanintersect = remove_rows_with_no_overlaps(intersect,'intersect_primary') # may want to wait, and make graphs for those with uces vs those without
			cleanintersect.rename(columns={'size_x':'size_{0}'.format(sfile)},inplace=True)
			cleanintersect.drop(columns=['size_y'],inplace=True)
			
			intersectstats = panda_describe_multiple_column(cleanintersect)
			
			allstats = pd.concat([intersectstats,primarystats,tertiarystats],axis=1)
			save_panda(allstats,'stats_{0}_intersections.txt'.format(primaryfile))
			
			# 6) generate graph for bin results
			cleanbin = drop_primary_zero_list(groupfeatures)
			formatbin = format_binned_data_for_graphing(cleanbin,bins)
			graph_binned_regions(formatbin,primaryfile,sfile)

if __name__ == "__main__":
	main()