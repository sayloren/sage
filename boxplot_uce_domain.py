"""
Script to make boxplots for regions with vs without uces

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
secondary - file with list of domain filenames (give the size of the regions with/without uces)
tertiary - genes (optional - if included will print the number of elements in this file per domain with/without uces)
quinary - mouse uces in original uce coordinates (optional - if included will add another box for regions that are only part of this subset)

Out:
pdf file with each of the domains in a seperate subplot, and all as the final most subplot

To Do:
make the gene boxplots not explicit, but an argument to supply

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
from scipy import stats
import math

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",required=True,type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeature",required=False,type=str,help="the tertiary features file")# Genes
	parser.add_argument("-q","--quinaryfeature",required=False,type=str,help="the quinary elements file - a subset of the primary features")# Mouse UCEs
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
	return pdcoordinates

# run primary, tertiary, quinary overlaps
def run_overlaps_for_ptq_against_s(secondary,pfile,tfile,qfile):
	pdprimary = count_overlap_df(secondary,pfile,'{0}'.format(pfile)) # make the pandas data sets for the count overlaps
	concat = pdprimary # rename to preserve
	if tfile: # if optional arguments, add to panda
		pdtertiary = count_overlap_df(secondary,tfile,'tertiary')
		concat = concat.merge(pdtertiary,how='inner',on=['chr','start','end','size'])
	if qfile: # if optional arguments, add to panda
		pdquinary = count_overlap_df(secondary,qfile,'{0}'.format(qfile))
		concat = concat.merge(pdquinary,how='inner',on=['chr','start','end','size'])
	return concat

# move elements without any overlaps
def remove_rows_with_no_overlaps(overlaps,column):
	return overlaps[overlaps[column]!=0]

# chunk data into number of graphs per page
def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i + n]

# format the with/without regions for graphing
def format_with_without_data_for_boxplot(pdfeatures,column,pfile,qfile):
	dropzero = (pdfeatures.loc[pdfeatures['intersect_{0}'.format(pfile)] > 0])
	sumdrop = dropzero[column].reset_index(drop=False)
	sumdrop['region'] = 'With UCEs'
	keepzero = (pdfeatures.loc[pdfeatures['intersect_{0}'.format(pfile)] < 1])
	sumkeep = keepzero[column].reset_index(drop=False)
	sumkeep['region'] = 'Without UCEs'
	tertiarysum = pd.concat([sumdrop,sumkeep])
	if qfile:
		keepquin = (pdfeatures.loc[pdfeatures['intersect_{0}'.format(qfile)] > 0])
		sumquin = keepquin[column].reset_index(drop=False)
		sumquin['region'] = 'With Mouse UCEs'
		tertiarysum = pd.concat([tertiarysum,sumquin])
	tertiarysum.drop(columns=['index'],inplace=True)
	#https://stackoverflow.com/questions/8362792/how-do-i-shift-the-decimal-place-in-python
	if column == 'size':
		tertiarysum['size'] /= 1000.
	return tertiarysum

# set the number of lines to darken on the boxplot
def set_num_lines_to_color(qfile):
	if qfile:
		numboxes,numlines = 3,9
	else:
		numboxes,numlines = 2,8
	return numboxes,numlines

# convert panda to bedtool
def panda_to_bedtool(panda):
	arArFeatures = panda.values.tolist()
	btFeatures = get_bedtools_features(arArFeatures)
	return btFeatures

# get standard deviation, from ruth's random region script
def getPopSD(arObservedOverlaps):
	floatLen = float(len(arObservedOverlaps))
	floatMean = float(sum(arObservedOverlaps))/len(arObservedOverlaps)
	dSumOfSquares = sum([((float(number) - floatMean) ** 2) for number in arObservedOverlaps])
	dVariance = float(dSumOfSquares) / floatLen
	return math.sqrt(dVariance)

# ks test from ruth's random region script
def KSTest(aOverlapBP):
	mean = float(sum(aOverlapBP)) / len(aOverlapBP)
	sd = getPopSD(aOverlapBP)
	rvNormMatched = stats.norm.rvs(loc=mean, scale=sd, size=len(aOverlapBP))
	npArOverlapBP = np.array(aOverlapBP)
	ksStat, KsPval = stats.ks_2samp(npArOverlapBP, rvNormMatched)
	if KsPval <= 0.05:
		strKSresult = "No"
	else:
		strKSresult = "Yes"
	return ksStat, KsPval, strKSresult

# run ks test for normal distribution and choose appropriate stats test
def run_appropriate_test(pdgroup,yvalue):
	withuces = pdgroup[yvalue].loc[pdgroup['region']=='With UCEs']
	withoutuces = pdgroup[yvalue].loc[pdgroup['region']=='Without UCEs']
	ksStat,KsPval,strKSresult = KSTest(withuces)
# 	if strKSresult == 'Yes':
# 		statcoef,statpval = stats.ttest_ind(withuces,withoutuces)# or ttest_rel()
# 		stattest = 'TT'
# 		formatpval = '{:.01e}'.format(statpval)
# 	else:
# 		statcoef,statpval = stats.mannwhitneyu(withuces,withoutuces)
# 		stattest = 'MW'
# 		formatpval = '{:.01e}'.format(statpval)
	statcoef,statpval = stats.mannwhitneyu(withuces,withoutuces)
	stattest = 'Mann Whiteny U p-value'
	formatpval = '{:.01e}'.format(statpval)
	return formatpval,stattest

# get the location where to add the p value annotation
def set_pval_label_location(pdgroup,yvalue):
# 	if yvalue == 'size':
# 		ylabelmax = pdgroup[yvalue].quantile(q=.99)+2
# 	else:
# 		ylabelmax = pdgroup[yvalue].quantile(q=.97)+2
	justelemenst = pdgroup.loc[pdgroup['region']=='With UCEs']
	excludeoutliers = justelemenst[np.abs(justelemenst[yvalue]-justelemenst[yvalue].mean())<=(3*justelemenst[yvalue].std())]
	ylabelmax = excludeoutliers[yvalue].max() + 2
	return ylabelmax

# darken the lines around the boxplot to black
def darkend_boxplot_lines(axes,numboxes,numlines,boxcolor):
	#https://stackoverflow.com/questions/36874697/how-to-edit-properties-of-whiskers-fliers-caps-etc-in-seaborn-boxplot
	for t,artist in enumerate(axes.artists):
		artist.set_edgecolor(boxcolor)
		for s in range(t*numboxes,t*numboxes+numlines):
			line = axes.lines[s]
			line.set_color(boxcolor)
			line.set_mfc(boxcolor)
			line.set_mec(boxcolor)

# tile the boxplots
def run_tiled_subplots_per_boxplot_dataset(pddata,yvalue,ylabeltext,names,filename,pfile,qfile):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(10,10))
# 	plt.rcParams['axes.formatter.limits'] = (-3, 3)
	sns.set_palette("Blues")
	datasetcounter = 0
	fig,ax_array = plt.subplots(3,2)
	intnum = len(names)
	numboxes,numlines = set_num_lines_to_color(qfile)
	for data_chunk,name_chunk in zip(chunks(pddata,6),chunks(names,6)):
		intPlotCounter = -1
		for i,ax_row in enumerate(ax_array):
			for j,axes in enumerate(ax_row):
				axes.cla()
				boxcolor = '#000000'
				intPlotCounter += 1
				if datasetcounter < len(names):
					pdgroup = format_with_without_data_for_boxplot(data_chunk[intPlotCounter],yvalue,pfile,qfile)
					sns.boxplot(data=pdgroup,x='region',y=yvalue,showfliers=False,ax=axes,linewidth=.75)
					axes.set_ylabel(ylabeltext,size=12)
					axes.set_xlabel('Domain Type',size=12)
					darkend_boxplot_lines(axes,numboxes,numlines,boxcolor)
					axes.set_title(name_chunk[intPlotCounter].split('.',1)[0],size=8)
					axes.set_xticklabels(axes.get_xticklabels(),fontsize=8)
					plt.setp(axes.xaxis.get_majorticklabels())#rotation=15
					formatpval,stattest = run_appropriate_test(pdgroup,yvalue)
					ylabelmax = set_pval_label_location(pdgroup,yvalue)
					axes.plot([0,0,1,1], [ylabelmax, ylabelmax+2, ylabelmax+2, ylabelmax], lw=.75, c=boxcolor)
					axes.text((0+1)*.5, ylabelmax+2,'{0}: {1}'.format(stattest,formatpval),ha='center',va='bottom',color=boxcolor,size=6,clip_on=False)
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
	pfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	if args.tertiaryfeature:
		tfile = args.tertiaryfeature
	else:
		tfile = None
	if args.quinaryfeature:
		qfile = args.quinaryfeature
	else:
		qfile = None

	# initiate collections
	lumpsecondary = []
	lumpsecondarycoords = []
	
	# process feature files
	for sfile in secondaryfiles:
		
		# get secondary features
		secondary = get_bedtools_features(sfile)
		
		# run count overlaps for primary, tertiary, quinary against the secondary regions
		concat = run_overlaps_for_ptq_against_s(secondary,pfile,tfile,qfile)
		
		# just the secondary coords
		coords = concat[['chr','start','end']]
		
		# add the coords to the lump later run a set for all the domains
		lumpsecondarycoords.append(coords)
		
		# add the data set to the lump sum to get the total for the end of the script
		lumpsecondary.append(concat)
	
	# concat lumped regions
	concatsecondary = pd.concat(lumpsecondary)
	concatsecondarycoords = pd.concat(lumpsecondarycoords)
	
	# convert lump secondary coordinates panda to bedtool
	btconcatsecondary = panda_to_bedtool(concatsecondarycoords)
	
	# run overlap with primary/tertiary/quinary again
	pdoverlapssecondary = run_overlaps_for_ptq_against_s(btconcatsecondary,pfile,tfile,qfile)
	
	# randomly select a subset of the total set based on how many sets of secondary regions input
	fracrandom = 1.0/len(secondaryfiles)
	pdrandomsecondary = pdoverlapssecondary.sample(frac=fracrandom)
	
	# add the concated all domains to the list to graph
	lumpsecondary.append(pdrandomsecondary)
	
	# add a descriptor to the concated domain dataset
	secondaryfiles.append('All Domains')
	
	# run the tile plot secondary sizes
	run_tiled_subplots_per_boxplot_dataset(lumpsecondary,'size','Size (kp)',secondaryfiles,'tiled_domain_sizes.pdf',pfile,qfile)
	
	if tfile:
		# run tile plot for tertiary counts
		run_tiled_subplots_per_boxplot_dataset(lumpsecondary,'intersect_tertiary','Frequency',secondaryfiles,'tiled_gene_number.pdf',pfile,qfile)

if __name__ == "__main__":
	main()