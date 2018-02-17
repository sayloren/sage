"""
Script to make point plots for the number of uces (or other feature) in each binned section of a domain

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

Question: only plot genes that are in regions that also have uces?
Question: plot the mouse uces against the total set?

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
import numpy as np

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",required=True,type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-b","--binnumber",type=int,default='10',help='number of bins to chunk the secondary files into, must be even number')
	return parser.parse_args()

# get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

# intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	return afile.intersect(bfile,c=True)

# convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	return pd.read_table(btfeature.fn,header=None)

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
	pdcoordinates.insert(len(pdcoordinates.columns),'id',range(0,0+len(pdcoordinates)))
	return pdcoordinates

# bin secondary regions
def make_window_with_secondary_files(sfile,bins):
	return pbt.BedTool().window_maker(b=sfile,n=bins,i="src")

# convert panda to bedtool
def convert_panda_to_bed_format(panda):
	arArFeatures = panda.values.tolist()
	return pbt.BedTool(arArFeatures)

# 2f) merge a and b files on 'id' column
def intersect_pandas_with_id(afile,bfile):
	return pd.merge(afile,bfile,on='id')

# 4e) groupby larger region
def group_df_by_secondary_regions(pdfeature):
	countcolumns = [col for col in pdfeature.columns if 'count' in col]
	outgroup = []
	for col in countcolumns:
		group = pd.DataFrame({'bincounts':pdfeature.groupby(['region_chr','region_start','region_end'])[col].apply(list)}).reset_index()
		group.columns = ['chr','start','end','bincounts_{0}'.format(col)]
		outgroup.append(group)
	return reduce(lambda x, y: pd.merge(x,y,on=['chr','start','end']),outgroup)

# remove the rows where there are no primary element overlaps with the secondary regions
def drop_primary_zero_list(pdfeatures,column):
	thresh = pdfeatures[~pdfeatures[column].apply(lambda row: all(item ==0 for item in row))]
	return thresh

# format the binned data frame for pointplot graphing
def format_binned_data_sum_for_graphing(pdfeatures,bins):
	selectcols = [col for col in pdfeatures.columns if 'bincounts' in col]
	listsum = []
	for group in selectcols:
		subset = pdfeatures[[group]]
		split = pd.DataFrame(subset[group].values.tolist())
		sum = split.sum(axis=0)
		listsum.append(sum)
	concatsum = pd.concat(listsum,axis=1)
	concatsum.columns = [selectcols]
	bincolumns = range(bins)
	format = pd.melt(concatsum)
	format.columns = ['filename','sumbin']
	format['bin'] = bincolumns * (format.shape[0]/len(bincolumns))
	return format

# fold data set sums
def fold_formated_binned_data_sum(pdfeatures,bins):
	halfbin = bins/2
	pdfeatures['invertsumbin'] = pdfeatures['sumbin'].iloc[::-1].reset_index(drop=True)
	pdfeatures['sumsums'] = pdfeatures['invertsumbin'] + pdfeatures['sumbin']
	headfeatures = pdfeatures.head(n=halfbin)
	dropfeatures = headfeatures[['filename','sumsums','bin']]
	dropfeatures.columns = ['filename','sumbin','bin']
	return dropfeatures

# chunk data into number of graphs per page
def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i + n]

# tile the point plots
def run_tiled_subplots_per_binned_dataset(pddata,names,filename):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(10,10))
	datasetcounter = 0
	fig,ax_array = plt.subplots(3,2)
	intnum = len(names)
	for data_chunk,name_chunk in zip(chunks(pddata,6),chunks(names,6)):
		intPlotCounter = -1
		for i,ax_row in enumerate(ax_array):
			for j,axes in enumerate(ax_row):
				axes.cla()
				intPlotCounter += 1
				if datasetcounter < len(names):
					pdgroup = data_chunk[intPlotCounter]
					sns.pointplot(data=pdgroup,x='bin',y='sumbin',color='#9ecae1',scale=1,ax=axes)
					axes.set_ylabel('Frequency',size=12)
					axes.set_xlabel('Bin Distance from Edge')
					axes.set_title(name_chunk[intPlotCounter].split('.',1)[0],size=8)
					for item in ([axes.xaxis.label] + axes.get_xticklabels()):
						item.set_fontsize(8)
					plt.setp(axes.xaxis.get_majorticklabels(),rotation=15)
					datasetcounter += 1
				else:
					axes.remove()
					pass
		plt.tight_layout()
		sns.despine()
		plt.savefig(pp, format='pdf')
	plt.clf()
	pp.close()

def main():
	args = get_args()
	
	# read in files from args
	bins = args.binnumber
	pfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	
	# initiate collection
	lumpsecondary = []
	
	# process feature files
	for sfile in secondaryfiles:
		
		# get secondary features
		secondary = get_bedtools_features(sfile)
		
		# make the pandas data sets for the count overlaps
		pdsecondary = count_overlap_df(secondary,pfile,'{0}'.format(pfile))
		
		# get just the coordinates and the id
		labelsecondary = pdsecondary[['chr','start','end','id']]
		
		# convert back to a bedtools to run the windows
		btsecondary = convert_panda_to_bed_format(labelsecondary)
		
		# run the domains through bedtools window to get the right # bins for each section
		windows = make_window_with_secondary_files(btsecondary,bins)
		
		# get the count of the number of uces in each window
		intersectprimary = intersect_bedfiles_c_true(windows,pfile)
		
		# convert to back to panda and label
		pdintersect = convert_bedtools_to_panda(intersectprimary)
		pdlabel = label_coordinate_columns(pdintersect)
		pdlabel.columns.values[3]='id'
		pdlabel.columns.values[4]='count_primary'
		pdlabel.rename(columns={'chr':'window_chr','start':'window_start','end':'window_end','size':'overlapsize'},inplace=True)
		
		# intersect back into domains
		intersectsecondary = intersect_pandas_with_id(pdlabel,labelsecondary)
		intersectsecondary.rename(columns={'chr':'region_chr','start':'region_start','end':'region_end'},inplace=True)
		intersectsecondary['region_size'] = intersectsecondary['region_end'].astype(int)-intersectsecondary['region_start'].astype(int)
		
		# group by the larger genomic regions; the whole domain
		groupsecondary = group_df_by_secondary_regions(intersectsecondary)
		
		# remove secondary features where there are no primary elements
		dropzero = drop_primary_zero_list(groupsecondary,'bincounts_count_primary')
		
		# shape the data frame
		formatbin = format_binned_data_sum_for_graphing(dropzero,bins)
		
		# fold the data frame by combining edges for the pointplot
		foldbin = fold_formated_binned_data_sum(formatbin,bins)
		
		# add the data set to the lump sum to get the total for the end of the script
		lumpsecondary.append(foldbin)
		
	# concat the lumped regions
	concatsecondary = pd.concat(lumpsecondary)
	
	# add the concated all domains to the list to graph
	lumpsecondary.append(concatsecondary)
	
	# add a descriptor to the concated domain dataset
	secondaryfiles.append('All Domains')
	
	# run tile plot for primaries binned
	run_tiled_subplots_per_binned_dataset(lumpsecondary,secondaryfiles,'tiled_binned_UCEs.pdf')







if __name__ == "__main__":
	main()