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
	sumdrop['region'] = 'Regions With UCEs'
	keepzero = (pdfeatures.loc[pdfeatures['intersect_{0}'.format(pfile)] < 1])
	sumkeep = keepzero[column].reset_index(drop=False)
	sumkeep['region'] = 'Regions Without UCEs'
	tertiarysum = pd.concat([sumdrop,sumkeep])
	if qfile:
		keepquin = (pdfeatures.loc[pdfeatures['intersect_{0}'.format(qfile)] > 0])
		sumquin = keepquin[column].reset_index(drop=False)
		sumquin['region'] = 'Regions With Mouse UCEs'
		tertiarysum = pd.concat([tertiarysum,sumquin])
	tertiarysum.drop(columns=['index'],inplace=True)
	return tertiarysum

# tile the boxplots
def run_tiled_subplots_per_boxplot_dataset(pddata,yvalue,ylabeltext,names,filename,pfile,qfile):
	sns.set_style('ticks')
	pp = PdfPages(filename)
	plt.figure(figsize=(10,10))
	plt.rcParams['axes.formatter.limits'] = (-3, 3)
	sns.set_palette("Blues")
	datasetcounter = 0
	fig,ax_array = plt.subplots(3,2)
	intnum = len(names)
	if qfile:
		numboxes,numlines = 3,9
	else:
		numboxes,numlines = 2,8
	for data_chunk,name_chunk in zip(chunks(pddata,6),chunks(names,6)):
		intPlotCounter = -1
		for i,ax_row in enumerate(ax_array):
			for j,axes in enumerate(ax_row):
				axes.cla()
				intPlotCounter += 1
				if datasetcounter < len(names):
					pdgroup = format_with_without_data_for_boxplot(data_chunk[intPlotCounter],yvalue,pfile,qfile)
					sns.boxplot(data=pdgroup,x='region',y=yvalue,showfliers=False,ax=axes,linewidth=.75)
					axes.set_ylabel(ylabeltext,size=12)
					axes.set_xlabel('')
# 					#https://stackoverflow.com/questions/36874697/how-to-edit-properties-of-whiskers-fliers-caps-etc-in-seaborn-boxplot
					for t,artist in enumerate(axes.artists):
						artist.set_edgecolor('#000000')
						for s in range(t*numboxes,t*numboxes+numlines):
							line = axes.lines[s]
							line.set_color('#000000')
							line.set_mfc('#000000')
							line.set_mec('#000000')
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

# def KSTest(aOverlapBP):
#     "Returns the KS test statistic and p value for rejecting the null hypothesis that aOverlapBP follows a normal distribution with mean and sd equal to those of aOverlapBP"
#     mean = float(sum(aOverlapBP)) / len(aOverlapBP)
#     if len(aOverlapBP) < 1000:
#         print 'Warning: number of iterations is < 1000; KS statistic may be unreliable'
#     sd = getPopSD(aOverlapBP)
#     rvNormMatched = stats.norm.rvs(loc=mean, scale=sd, size=len(aOverlapBP))
#     npArOverlapBP = np.array(aOverlapBP)
#     ksStat, KsPval = stats.ks_2samp(npArOverlapBP, rvNormMatched)
#     if KsPval <= 0.05:
#         strKSresult = "No"
#         print 'KS statistic is significant: attention needed'
#     else:
#         strKSresult = "Yes"
#         print 'KS statistic not significant: random overlaps appear normally distributed'
#     return ksStat, KsPval, strKSresult

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

	# initiate collection
	lumpsecondary = []
	
	# process feature files
	for sfile in secondaryfiles:
		
		# get secondary features
		secondary = get_bedtools_features(sfile)
		
		# make the pandas data sets for the count overlaps
		pdprimary = count_overlap_df(secondary,pfile,'{0}'.format(pfile))
		
		# rename to preserve
		concat = pdprimary
		
		# if optional arguments, add to panda
		if tfile:
			pdtertiary = count_overlap_df(secondary,tfile,'tertiary')
			concat = concat.merge(pdtertiary,how='inner',on=['chr','start','end','size'])
		if qfile:
			pdquinary = count_overlap_df(secondary,qfile,'{0}'.format(qfile))
			concat = concat.merge(pdquinary,how='inner',on=['chr','start','end','size'])
		
		# add the data set to the lump sum to get the total for the end of the script
		lumpsecondary.append(concat)
	
	# concat the lumped regions
	concatsecondary = pd.concat(lumpsecondary)
	
	# add the concated all domains to the list to graph
	lumpsecondary.append(concatsecondary)
	
	# add a descriptor to the concated domain dataset
	secondaryfiles.append('All Domains')
	
	# run the tile plot secondary sizes
	run_tiled_subplots_per_boxplot_dataset(lumpsecondary,'size','Size (bp)',secondaryfiles,'tiled_domain_sizes.pdf',pfile,qfile)
	
	# run tile plot for tertiary counts
	run_tiled_subplots_per_boxplot_dataset(lumpsecondary,'intersect_tertiary','Frequency',secondaryfiles,'tiled_gene_number.pdf',pfile,qfile)

if __name__ == "__main__":
	main()