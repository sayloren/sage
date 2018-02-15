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
import numpy as np

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",required=True,type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeatures",type=argparse.FileType('rU'),help="a file with a list of file names with the tertiary features to query")# Genes
	parser.add_argument("-q","--quinaryfeatures",type=str,help="the quinary elements file - a subset of the primary features")# Mouse UCEs
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

# 4d) 7 cols coord labels
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
	primarylabel.columns.values[4]='count_primary'
	primarylabel.rename(columns={'chr':'achr','start':'astart','end':'aend','size':'overlapsize'},inplace=True)
	sfile.rename(columns={'chr':'bchr','start':'bstart','end':'bend'},inplace=True)
	intersectsecondary = intersect_pandas_with_id(primarylabel,sfile)
	outintersect = intersectsecondary[['achr','astart','aend','count_primary','bchr','bstart','bend','overlapsize','id']]
	return outintersect,windows

# 3) Query: What other features characterize domains with UCEs; in a domain
def additional_features_overlap(secondary,scoord,tfile):
	tertiary = get_bedtools_features(tfile)
	intersect = intersect_bedfiles_c_true(secondary,tertiary)
	pdfeature = convert_bedtools_to_panda(intersect)
	tcoord = label_coordinate_columns(pdfeature)
	tcoord.columns.values[3]='intersect_{0}'.format(tfile)
	concatintersect = pd.merge(scoord,tcoord,how='inner',right_on=['chr','start','end'],left_on=['bchr','bstart','bend'])
	return tertiary,concatintersect

# 4) Query: What other features characterize domains with UCEs; where in domain
def additional_features_locate(secondary,tertiary,tfile,scoord,windows):
	intersecttertiary = intersect_bedfiles_c_true(windows,tertiary)
	tertiarypd = convert_bedtools_to_panda(intersecttertiary)
	tertiarylabel = label_coordinate_columns(tertiarypd)
	tertiarylabel.columns.values[3]='id'
	tertiarylabel.columns.values[4]='count_{0}'.format(tfile)
	tertiarylabel.rename(columns={'chr':'cchr','start':'cstart','end':'cend','size':'overlapsize'},inplace=True)
	concatintersect = pd.merge(scoord,tertiarylabel,left_on=['achr','astart','aend','id','overlapsize'],right_on=['cchr','cstart','cend','id','overlapsize'])#,how='inner'
	concatintersect.drop(['cchr','cstart','cend'],axis=1,inplace=True)
	pdgroup = group_df_by_secondary_regions(concatintersect)
	return pdgroup

# remove the rows where there are no primary element overlaps with the secondary regions
def drop_primary_zero_list(pdfeatures,column):
	thresh = pdfeatures[~pdfeatures[column].apply(lambda row: all(item ==0 for item in row))]
	return thresh

# keep only those rows where there are no primary element overlaps with the secondary regions
def keep_primary_zero_list(pdfeatures,column):
	thresh = pdfeatures[pdfeatures[column].apply(lambda row: all(item ==0 for item in row))]
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

# format the binned data frame for box plot graphing
def format_binned_dataframe_for_boxplot(pdfeatures,bins):
	selectcols = [col for col in pdfeatures.columns if 'bincounts' in col]
	listdf = []
	for group in selectcols:
		subset = pdfeatures[[group]]
		split = pd.DataFrame(subset[group].values.tolist())
		split['filename'] = group
		listdf.append(split)
	concatdf = pd.concat(listdf,axis=0)
	return concatdf

# format the with/without regions for graphing
def format_with_without_data_for_boxplot(pdfeatures,column,quinaryfiles):
	dropzero = (pdfeatures.loc[pdfeatures['intersect_primary'] > 0])
	sumdrop = dropzero[column].reset_index(drop=False)
	sumdrop['primary'] = 'Regions With UCEs'
	keepzero = (pdfeatures.loc[pdfeatures['intersect_primary'] < 1])
	sumkeep = keepzero[column].reset_index(drop=False)
	sumkeep['primary'] = 'Regions Without UCEs'
	keepquin = (pdfeatures.loc[pdfeatures['intersect_{0}'.format(quinaryfiles)] > 0])
	sumquin = keepquin[column].reset_index(drop=False)
	sumquin['primary'] = 'Regions With Mouse UCEs'
	tertiarysum = pd.concat([sumdrop,sumkeep,sumquin])
	tertiarysum.drop(columns=['index'],inplace=True)
	return tertiarysum

# fold data set sums
def fold_formated_binned_data_sum(pdfeatures,bins):
	halfbin = bins/2
	pdfeatures['invertsumbin'] = pdfeatures['sumbin'].iloc[::-1].reset_index(drop=True)
	pdfeatures['sumsums'] = pdfeatures['invertsumbin'] + pdfeatures['sumbin']
	headfeatures = pdfeatures.head(n=halfbin)
	dropfeatures = headfeatures[['filename','sumsums','bin']]
	dropfeatures.columns = ['filename','sumbin','bin']
	return dropfeatures

# fold data set 
def fold_formated_binned_data_df(pdfeatures,bins):
	filefeatures = pdfeatures[['filename']]
	pdfeatures.drop('filename',axis=1,inplace=True)
	halfbin = bins/2
	reversefeatures = pdfeatures.iloc[:, ::-1]
	reversefeatures.columns = range(bins)
	sumfeatures = pdfeatures.add(reversefeatures, fill_value=0)
	halffeatures = sumfeatures.iloc[:,:halfbin]
	catfeatures = pd.concat([halffeatures,filefeatures],axis=1)
	return catfeatures

# line graph binned data with no tertiary features
def graph_binned_regions_no_tertiary(pdfeatures,filename):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(14,7))
	
	unique = len(pdfeatures['filename'].unique())
	sns.set_palette("Blues",n_colors=unique)
	
	gs = gridspec.GridSpec(1,1)
	gs.update(hspace=.8)
	
	ax0 = plt.subplot(gs[0,0])
	sns.pointplot(data=pdfeatures,x='bin',y='sumbin',color='#9ecae1',scale=3)
	ax0.set_ylabel('Frequency')
	ax0.set_xlabel('Bin Distance from Edge')
	for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] + ax0.get_xticklabels() + ax0.get_yticklabels()):
		item.set_fontsize(22)
	
	sns.despine()
	pp.savefig()
	pp.close()

# normalize data
def normalize_by_total_count(bins,df,filename,length):
	halfbin = bins/2
	subsetprimaryfiles = df[df['filename']==filename]
	subsetprimarydata = subsetprimaryfiles.iloc[:,:halfbin]
	normalizeprimary = subsetprimarydata.apply(lambda x: x/length)
	normalizeprimary['filename'] = filename
	return normalizeprimary

# graph box plots for secondary and tertiary data - will have to put in the y axis label as arg
def graph_boxplot_region_size(pdfeatures,filename,yvalue,ylabeltext):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(14,7))
	plt.rcParams['axes.formatter.limits'] = (-3, 3)
	
	sns.set_palette("Blues")
	
	gs = gridspec.GridSpec(1,1)
	gs.update(hspace=.8)
	
	ax0 = plt.subplot(gs[0,0])
	sns.boxplot(data=pdfeatures,x='primary',y=yvalue,showfliers=False)
	ax0.set_ylabel(ylabeltext)
	ax0.set_xlabel('')
	for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] + ax0.get_xticklabels() + ax0.get_yticklabels()):
		item.set_fontsize(22)
	sns.despine()
	pp.savefig()
	pp.close()

# graph box plots for secondary and tertiary data over binned regions
def graph_boxplot_binned_regions(pdfeatures,filename):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(14,7))
	
	sns.set_palette("Blues")
	
	gs = gridspec.GridSpec(1,1)
	gs.update(hspace=.8)
	
	ax0 = plt.subplot(gs[0,0])
	sns.boxplot(data=pdfeatures,x='variable',y='value',showfliers=False,hue='filename')
	ax0.set_ylabel('Fraction of Total Element Count')
	ax0.set_xlabel('Bin Distance from Edge')
	for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] + ax0.get_xticklabels() + ax0.get_yticklabels()):
		item.set_fontsize(22)
	
	sns.despine()
	pp.savefig()
	pp.close()

# get the number of quinary features
def number_quinary_features(quinaryfiles):
	quinary = get_bedtools_features(quinaryfiles)
	pdquinary = convert_bedtools_to_panda(quinary)
	lengthquinary = len(pdquinary)
	print 'there are {0} elements in {1}'.format(lengthquinary,quinaryfiles)
	return quinary

# intersect quinary features to get the subset of primary featutes
def intersect_quinary_features(quinary,scoord,quinaryfiles,sfile):
	tempsecondary,qcoord = overlaping_features(quinary,sfile)
	qcoord.rename(columns={'intersect_primary':'intersect_{0}'.format(quinaryfiles)},inplace=True)
	qcoord.drop(labels=['id'],axis=1,inplace=True)
	concatintersect = scoord.merge(qcoord,how='left',on=['chr','start','end','size'])
	return qcoord,concatintersect

# print number of features from secondary intersections with primary and quinary
def print_number_features(scoord,qcoord,sfile,quinaryfiles):
	sumuce = scoord['intersect_primary'].sum()
	summouse = qcoord['intersect_{0}'.format(quinaryfiles)].sum()
	print '{0} total uces, {1} mouse uces that are overlapped in secondary file {2}'.format(sumuce,summouse,sfile)

# rename columns in stats file for clarity
def rename_stats_columns(alltertiarystats,primaryfile,sfile,tfile,quinaryfiles):
	alltertiarystats.rename(columns={'intersect_primary':'intersect_{0}_{1}'.format(primaryfile,sfile)},inplace=True)
	alltertiarystats.rename(columns={'intersect_{0}'.format(tfile):'intersect_{0}_{1}'.format(tfile,sfile)},inplace=True)
	alltertiarystats.rename(columns={'intersect_{0}'.format(quinaryfiles):'intersect_{0}_{1}'.format(quinaryfiles,sfile)},inplace=True)
	alltertiarystats.rename(columns={'size'.format(tfile):'size_{0}'.format(sfile)},inplace=True)
	alltertiarystats.rename(columns={'size_primary'.format(tfile):'size_{0}'.format(primaryfile)},inplace=True)
	return alltertiarystats

# tile the plots for all features
def tile_all_collected_features(graphone,graphtwo,graphthree,graphfour,filename):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(12,7))
	
	sns.set_palette("Blues")
	
	gs = gridspec.GridSpec(2,2)
	gs.update(hspace=.5)
	
	ax0 = plt.subplot(gs[0,0])
	sns.pointplot(data=graphone,x='bin',y='sumbin',color='#9ecae1',scale=3)
	ax0.set_ylabel('Frequency')
	ax0.set_xlabel('Bin Distance from Edge')
	for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] + ax0.get_xticklabels() + ax0.get_yticklabels()):
		item.set_fontsize(22)

	ax1 = plt.subplot(gs[0,1])
	sns.pointplot(data=graphtwo,x='bin',y='sumbin',color='#9ecae1',scale=3)
	ax1.set_ylabel('Frequency')
	ax1.set_xlabel('Bin Distance from Edge')
	for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
		item.set_fontsize(22)

	ax2 = plt.subplot(gs[1,0])
	sns.boxplot(data=graphthree,x='primary',y='size',showfliers=False)
	ax2.set_ylabel('Size (bp)')
	ax2.set_xlabel('')
	for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_yticklabels()):
		item.set_fontsize(22)
	plt.setp(ax2.xaxis.get_majorticklabels(),rotation=15)
	
	ax3 = plt.subplot(gs[1,1])
	sns.boxplot(data=graphfour,x='primary',y='intersect_tertiary',showfliers=False)
	ax3.set_ylabel('Frequency')
	ax3.set_xlabel('')
	for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] + ax3.get_yticklabels()):
		item.set_fontsize(22)
	plt.setp(ax3.xaxis.get_majorticklabels(),rotation=15)
	sns.despine()
	pp.savefig()
	pp.close()

# run the graphs for the lumped all datasets
def lumped_all_data_set_graphs(allstats,allsecsize,allbinprime,allbintert,allsumtert,pfile,qfiles):
	# concat all the primary elements in the binned secondary regions
	catprimarybin = pd.concat(allbinprime)
	# graph the primary binned elements for all secondary features
# 	graph_binned_regions_no_tertiary(catprimarybin,'bincounts_uces_all.pdf')
	# concat all the tertiary elements in the binned secondary regions
	cattertiarybin = pd.concat(allbintert)
	# graph the binned tertiary elements for all secondary features
# 	graph_binned_regions_no_tertiary(cattertiarybin,'bincounts_genes_all.pdf')
	# concat all the tertiary sums
	cattertiarysum = pd.concat(allsumtert)
	totaltert = format_with_without_data_for_boxplot(cattertiarysum,'intersect_tertiary',qfiles)
	# graph the tertiary sums in box plot
# 	graph_boxplot_region_size(totaltert,'genecount_all.pdf','intersect_tertiary','Frequency')
	# make the stats data frame for all the secondary features
	concatsecondary = pd.concat(allstats)
	concatsecondary.rename(columns={'size':'size_all_secondary'},inplace=True)
	catsecondarystats = panda_describe_multiple_column(concatsecondary)
	# format data for boxplot graphs of secondary size
	concatsecondaryfull = pd.concat(allsecsize)
	totalsec = format_with_without_data_for_boxplot(concatsecondaryfull,'size',qfiles)
	# graph counts for all secondary feature sizes
# 	graph_boxplot_region_size(totalsec,'domainsize_all.pdf','size','Size (bp)')
	# save stats to file
	save_panda(catsecondarystats,'stats_{0}_intersection.txt'.format(pfile))
	# tile all the concated 'all' graphs
	tile_all_collected_features(catprimarybin,cattertiarybin,totalsec,totaltert,'tiled_all_graphs.pdf')

# tile the boxplots
def run_tiled_subplots_per_boxplot_dataset(pddata,yvalue,ylabeltext,names,filename,quinaryfiles):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(10,10))
	plt.rcParams['axes.formatter.limits'] = (-3, 3)
	sns.set_palette("Blues")
	plt.tight_layout()
	sns.despine()

# 	if intnum % 2 != 0:
# 		fig.delaxes(ax_array[intnum / 2, 1])

	datasetcounter = 0
	fig,ax_array = plt.subplots(3,2)
	intnum = len(names)
	for data_chunk,name_chunk in zip(chunks(pddata,6),chunks(names,6)):
		intPlotCounter = -1
		for i,ax_row in enumerate(ax_array):
			for j,axes in enumerate(ax_row):
				intPlotCounter += 1
				if datasetcounter < len(names):
					pdgroup = format_with_without_data_for_boxplot(data_chunk[intPlotCounter],yvalue,quinaryfiles)
					sns.boxplot(data=pdgroup,x='primary',y=yvalue,showfliers=False,ax=axes)
					axes.set_ylabel(ylabeltext,size=12)
					axes.set_xlabel('')
					axes.set_title(name_chunk[intPlotCounter],size=8)
					for item in ([axes.xaxis.label] + axes.get_xticklabels()):
						item.set_fontsize(8)
					plt.setp(axes.xaxis.get_majorticklabels(),rotation=15)
					datasetcounter += 1
				else:
					pass
		plt.savefig(pp, format='pdf')
	plt.clf()
	pp.close()

# chunk data into number of graphs per page
def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i + n]

# tile the point plots
def run_tiled_subplots_per_binned_dataset(pddata,names,filename):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(10,10))
	ncols = 2
	intnum = len(names)
	if intnum % 2 ==0:
		nrows = intnum/2
		fig,ax_array = plt.subplots(nrows,ncols)
	elif intnum % 2 != 0:
		nrows = (intnum/2) + 1
		fig,ax_array = plt.subplots(nrows,ncols)
	intPlotCounter = -1
	for i,ax_row in enumerate(ax_array):
		for j,axes in enumerate(ax_row):
			intPlotCounter += 1
			if intPlotCounter < len(pddata):
				pdgroup = pddata[intPlotCounter]
				sns.pointplot(data=pdgroup,x='bin',y='sumbin',color='#9ecae1',scale=1,ax=axes)
				axes.set_ylabel('Frequency',size=12)
				axes.set_xlabel('Bin Distance from Edge')
				axes.set_title(names[intPlotCounter],size=8)
			else:
				pass
	if intnum % 2 != 0:
		fig.delaxes(ax_array[intnum / 2, 1])
	plt.tight_layout()
	sns.despine()
	pp.savefig()
	plt.clf()
	pp.close()

def main():
	args = get_args()
	
	# read in files from args
	primaryfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tertiaryfiles = [line.strip() for line in args.tertiaryfeatures]
	quinaryfiles = args.quinaryfeatures # have to integrate to graph subset of uces on boxplot
	bins = args.binnumber
	genomefile = args.genomefile
	
	# make list of all files to iterate through
	allfiles = secondaryfiles + tertiaryfiles
	allfiles.insert(0,primaryfile)
	
	# get the primary features
	primary = get_bedtools_features(primaryfile)
	
	# get the number of primary features
	pdprimary = convert_bedtools_to_panda(primary)
	lengthprimary = len(pdprimary)
	print 'there are {0} elements in {1}'.format(lengthprimary,primaryfile)
	
	# initiate collections
	lumpstats = []
	lumpsecondary = []
	lumpprimarybin = []
	lumptertiarybin = []
	lumptertiary = []
	
	lumptertiary = []
	
	# process feature files
	for sfile in secondaryfiles:
		# 1) Query: How many UCEs are there in a domain
		secondary,scoord = overlaping_features(primary,sfile)
		
		# get the number of secondary features with no primary overlaps
		nooverlaps = count_number_with_zero_overlaps(scoord,'intersect_primary')
		print '{0} instances of no overlaps of primary element on {1}'.format(nooverlaps,sfile)
		
		# print number of quinary features, return quinary features in panda
		quinary = number_quinary_features(quinaryfiles)
		
		# make quinary intersections
		qcoord,concatintersect = intersect_quinary_features(quinary,scoord,quinaryfiles,sfile)
		
		# get the number of uces in the secondary file
		print_number_features(scoord,qcoord,sfile,quinaryfiles)
		
		# remove all secondary regions with no primary overlaps
		cleanintersect = remove_rows_with_no_overlaps(concatintersect,'intersect_primary')
		
		# 2a) Query: Where in the domain are the UCEs
		binfeatures,windows = locateing_features(primary,scoord,bins)

		# group bin counts by secondary features they came from
		groupfeatures = group_df_by_secondary_regions(binfeatures)
		
		# remove secondary features where there are no primary elements
		dropzero = drop_primary_zero_list(groupfeatures,'bincounts_count_primary')
		
		# shape the data frame
		formatbin = format_binned_data_sum_for_graphing(dropzero,bins)
		
		# fold the data frame by combining edges for the pointplot
		foldbin = fold_formated_binned_data_sum(formatbin,bins)
		
		lumpprimarybin.append(foldbin)
		
		# generate graph for bin results
# 		graph_binned_regions_no_tertiary(foldbin,'bincounts_uces_{0}_{1}.pdf'.format(primaryfile,sfile))
		
		# format data for boxplot graphs of secondary size
		secondarysum = format_with_without_data_for_boxplot(concatintersect,'size',quinaryfiles)
		
		# graph size for secondary regions
# 		graph_boxplot_region_size(secondarysum,'domainsize_{0}.pdf'.format(sfile),'size','Size (bp)')
		
		# add the primary x secondary intersections data to the list for lump secondary stats
		lumpstats.append(cleanintersect)
		lumpsecondary.append(concatintersect)
		
		# for contrast with tertiary features
		for tfile in tertiaryfiles:
			# 3) Query: What other features characterize domains with UCEs; in a domain
			tertiary,intersect = additional_features_overlap(secondary,scoord,tfile)
			
			# concat quinary data set
			concatintersect = pd.merge(intersect,qcoord,how='left',right_on=['chr','start','end'],left_on=['bchr','bstart','bend'])
			concatintersect.drop(labels=['chr_x','size_x','start_x','end_x','chr_y','size_y','start_y','end_y'],axis=1,inplace=True)
			
			# get the number of tertiary feature
			pdtertiary = convert_bedtools_to_panda(tertiary)
			lengthtertiary = len(pdtertiary)
			
			# 4) Query: What other features characterize domains with UCEs; where in domain
			groupfeatures = additional_features_locate(secondary,tertiary,tfile,binfeatures,windows)
			
			# remove secondary features where there are no primary elements
			dropzero = drop_primary_zero_list(groupfeatures,'bincounts_count_primary')
			
			# format the data frame with the filename as a column
			formatdf = format_binned_dataframe_for_boxplot(dropzero,bins)
			formatbin = format_binned_data_sum_for_graphing(dropzero,bins)
			
			# fold the data frame by combining edges for the boxplot
			folddf = fold_formated_binned_data_df(formatdf,bins)
			
			tertiarybin = formatbin[formatbin['filename']=='bincounts_count_{0}'.format(tfile)].reset_index(drop=True)
			
			foldbin = fold_formated_binned_data_sum(tertiarybin,bins)
			
			lumptertiarybin.append(foldbin)
			
			# dot plot with gene size for uce-domain intersected regions
# 			graph_binned_regions_no_tertiary(foldbin,'bincounts_genes_{0}_{1}_{2}.pdf'.format(primaryfile,sfile,tfile))
			
			# normalize the values by total counts per each data set
			normalizedprimary = normalize_by_total_count(bins,folddf,'bincounts_count_primary',lengthprimary)
			normalizedtertiary = normalize_by_total_count(bins,folddf,'bincounts_count_{0}'.format(tfile),lengthtertiary)
			catnormalized = pd.concat([normalizedprimary,normalizedtertiary])
			
			meltnormalized = pd.melt(catnormalized,id_vars=['filename'])
			
			# graph normalized by total count elements - non-functional!
# 			graph_boxplot_binned_regions(meltnormalized,'bincount_dist_{0}_{1}.pdf'.format(primaryfile,sfile))
			
			# 5) generate stats results
			primarystats = panda_describe_single_column(primary,'size_primary')
			tertiarystats = panda_describe_single_column(tertiary,'size_{0}'.format(tfile))
			
			# get the number of features with no primary overlaps
			nooverlaps = count_number_with_zero_overlaps(intersect,'intersect_primary')
			print '{0} instances of no overlaps of primary element on {1}'.format(nooverlaps,sfile)
			
			# remove all secondary regions with no primary overlaps
			cleanintersect = remove_rows_with_no_overlaps(concatintersect,'intersect_primary')
			
			# make stats file
			intersectstats = panda_describe_multiple_column(cleanintersect)
			alltertiarystats = pd.concat([intersectstats,primarystats,tertiarystats],axis=1)
			
			# make sure names are clear
			alltertiarystats = rename_stats_columns(alltertiarystats,primaryfile,sfile,tfile,quinaryfiles)
			
			# save stats to file
			save_panda(alltertiarystats,'stats_{0}_intersection.txt'.format(primaryfile))
			
			# format data for boxplot graphs of tertiary feature counts
			tertiarysum = format_with_without_data_for_boxplot(concatintersect,'intersect_{0}'.format(tfile),quinaryfiles)
			
			tertiarysum.rename(columns={'intersect_{0}'.format(tfile):'intersect_tertiary'},inplace=True)
			concatintersect.rename(columns={'intersect_{0}'.format(tfile):'intersect_tertiary'},inplace=True)
			
# 			lumptertiary.append(tertiarysum)
			lumptertiary.append(concatintersect)
			
			# graph counts for tertiary features
# 			graph_boxplot_region_size(tertiarysum,'genecount_{0}_{1}.pdf'.format(sfile,tfile),'intersect_tertiary','Frequency')
	
	# run the tile plot secondary sizes
	run_tiled_subplots_per_boxplot_dataset(lumpsecondary,'size','Size (bp)',secondaryfiles,'tiled_domain_sizes.pdf',quinaryfiles)
	
	# run tile plot for tertiary counts
	run_tiled_subplots_per_boxplot_dataset(lumptertiary,'intersect_tertiary','Frequency',secondaryfiles,'tiled_gene_number.pdf',quinaryfiles)
	
	# run tile plot for primaries binned
	run_tiled_subplots_per_binned_dataset(lumpprimarybin,secondaryfiles,'tiled_binned_UCEs.pdf')
	
	# run tile plot for tertiary binned
	run_tiled_subplots_per_binned_dataset(lumptertiarybin,secondaryfiles,'tiled_binned_genes.pdf')
	
	# run the lump all graphs
	lumped_all_data_set_graphs(lumpstats,lumpsecondary,lumpprimarybin,lumptertiarybin,lumptertiary,primaryfile,quinaryfiles)

if __name__ == "__main__":
	main()