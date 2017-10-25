"""
Script to make karyogram pictures, where annotations can be either dense or spread out on each chromosome

Wren Saylor
August 21 2017
"""

import argparse
import pybedtools as pbt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import numpy as np
from collections import defaultdict

def get_args():
	parser = argparse.ArgumentParser(description="plot a karyogram for some bedfile data")
	parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files to annotate')
	parser.add_argument("-k","--karyogramdata", type=str, default="karyotype_hg19.txt")
	parser.add_argument("-c","--colorcolumn",type=int,default="3",help="column to use for annotations, must be the same for all annotation files, remember this is 0 based, so subtract 1")
	parser.add_argument("-a","--annotationcolorcolumn",type=int,default="4",help="column to use chromosome banding colors from the karyogram file, remember this is 0 based, so subtract 1")
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-s","--spread",action='store_true', help='add if you want each chromosome on a different page with the full spread of annotations')
	return parser.parse_args()

def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

def save_bedtool_as_bedfile(btObject,strFilename):
	btObject.saveas(strFilename)

def convert_bedtool_to_panda(btobject):
	save_bedtool_as_bedfile(btobject,'temp.bed')
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

def get_data(btobject):
	btfeatures = get_bedtools_features(btobject)
	pdfeatures = convert_bedtool_to_panda(btfeatures)
	return pdfeatures

def select_chromosomes(genome):
	chromosomeSet = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	keepSizes = genome[genome[0].isin(chromosomeSet)]
	keepSizes.columns = ['chromosome','size']
	sorterIndex = dict(zip(chromosomeSet,range(len(chromosomeSet))))
	keepSizes['chromosome_Rank'] = keepSizes['chromosome'].map(sorterIndex)
	keepSizes.sort_values([('chromosome_Rank')],inplace=True)
	keepSizes.drop('chromosome_Rank',1, inplace = True)
	keepSizes.reset_index(inplace=True,drop=True)
	return keepSizes

def unique_annotation_element_dictionary(annotationFiles,column):
	collectUnique = []
	for file in annotationFiles:
		pdfeatures = get_data(file)
		unique = pdfeatures[column].unique()
		uniqueList = unique.tolist()
		collectUnique.append(uniqueList)
	flatUnique = [item for sublist in collectUnique for item in sublist]
	colorPaletteAnnotation = sns.color_palette(palette="cubehelix",n_colors=len(flatUnique),desat=.9)
	colorDictAnnotation = dict(zip(flatUnique,colorPaletteAnnotation))
	return colorDictAnnotation

def set_plot_params_dense(annotationFiles,regions,column,annotationcolumn,sortedChrSect,length_first,ax,colorDictKaryogram,colorDictAnnotation):
	num_annotation_files = len(annotationFiles)
	annotation_space = .1
	DIM = 1.0
	order=1+(num_annotation_files*annotation_space*2)
	ax.set_xlim([0.0, DIM * (1.75)+(num_annotation_files*annotation_space)])
	ax.set_ylim([0.0, DIM])
	for index,chromosome in sortedChrSect.iterrows():
		chr_size = chromosome['size']
		first_size = length_first['size']
		x_start=order * DIM * 0.1
		x_end=x_start + (DIM * 0.06)
		y_start=DIM * 0.9 * (chr_size/first_size) 
		y_end=DIM * 0.1
		# rounded edges
		radius=(x_end-x_start)/2.0
		center_x=x_start+(x_end-x_start)/2.0
		subset_chr = regions.loc[regions[0]==chromosome['chromosome']]
		end_coord = subset_chr.loc[subset_chr[2].idxmax()]
		end_color = colorDictKaryogram[end_coord[annotationcolumn]]
		start_coord = subset_chr.loc[subset_chr[1].idxmin()]
		start_color = colorDictKaryogram[start_coord[annotationcolumn]]
		w1 = Wedge((center_x,y_start),radius,-90.0,270.0,facecolor=start_color,edgecolor='black',linewidth=1.5)
		w2 = Wedge((center_x,y_end),radius,-270.0,90.0,facecolor=end_color,edgecolor='black',linewidth=1.5)
		ax.add_patch(w1)
		ax.add_patch(w2)
		ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
		ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')
		# iterate through annotations for each chromosome
		for index,band in subset_chr.iterrows():
			current_height_band = band[2]-band[1]
			current_height_band_scaled=((y_end-y_start)/chr_size)*current_height_band
			band_scaled=y_start+((y_end-y_start)/chr_size)*band[1]
			band_color = colorDictKaryogram[band[annotationcolumn]]
			r=Rectangle((x_start,band_scaled),x_end-x_start,current_height_band_scaled,color=band_color)
			ax.add_patch(r)
		annotation_coord = (DIM*.015)
		for file in annotationFiles:
			pdfeatures = get_data(file)
			pd_chr = pdfeatures.loc[pdfeatures[0]==chromosome['chromosome']]
			for index,point in pd_chr.iterrows():
				current_height_point = point[2]-point[1]
				current_height_point_scaled=((y_end-y_start)/chr_size)*current_height_point
				point_scaled=y_start+((y_end-y_start)/chr_size)*point[1]
				point_color = colorDictAnnotation[point[column]]
				r=Rectangle(((x_end+annotation_coord),point_scaled),.01,current_height_point_scaled,color=point_color)
				ax.add_patch(r)
			annotation_coord = annotation_coord+(DIM*.015)
		ax.text(center_x-.045, y_end - (DIM * 0.07), chromosome['chromosome'])
		order=order+1+(num_annotation_files*annotation_space*2)
	collectPatches = []
	for key,value in colorDictAnnotation.iteritems():
		patch = patches.Patch(color=value, label=key)
		collectPatches.append(patch)
	plt.legend(handles=collectPatches)
	ax.set_yticklabels([])
	ax.set_xticklabels([])

def set_plot_params_spread(annotationFiles,column,annotationcolumn,chromosome,length,chromosome_label,ax,colorDictKaryogram,colorDictAnnotation):
	num_annotation_files = len(annotationFiles)
	annotation_space = .1
	DIM = 1.0
	ax.set_xlim([0.0, DIM * (1.75)+(num_annotation_files*annotation_space)])
	ax.set_ylim([0.0, DIM])
	chr_size = length['size']
	x_start=(num_annotation_files*annotation_space*2) * DIM * 0.1
	x_end=x_start + (DIM * 0.06)
	y_start=DIM * 0.9
	y_end=DIM * 0.1
	# rounded edges
	radius=(x_end-x_start)/2.0
	center_x=x_start+(x_end-x_start)/2.0
	end_coord=chromosome.loc[chromosome[2].idxmax()]
	end_color=colorDictKaryogram[end_coord[annotationcolumn]]
	start_coord=chromosome.loc[chromosome[1].idxmin()]
	start_color=colorDictKaryogram[start_coord[annotationcolumn]]
	w1=Wedge((center_x,y_start),radius,-90.0,270.0,facecolor=start_color,edgecolor='black',linewidth=1.5)#-.006
	w2=Wedge((center_x,y_end),radius,-270.0,90.0,facecolor=end_color,edgecolor='black',linewidth=1.5)#+.006
	ax.add_patch(w1)
	ax.add_patch(w2)
	ax.plot([x_start,x_start],[y_start,y_end],ls='-',color='black')
	ax.plot([x_end,x_end],[y_start,y_end],ls='-',color='black')
	# iterate through annotations for each chromosome
	for index,band in chromosome.iterrows():
		current_height_band=band[2]-band[1]
		current_height_band_scaled=((y_end-y_start)/chr_size)*current_height_band
		band_scaled=y_start+((y_end-y_start)/chr_size)*band[1]
		band_color=colorDictKaryogram[band[annotationcolumn]]
		r=Rectangle((x_start,band_scaled),x_end-x_start,current_height_band_scaled,color=band_color)
		ax.add_patch(r)
	annotation_coord=(DIM*.015)
	for file in annotationFiles:
		pdfeatures=get_data(file)
		pd_chr=pdfeatures.loc[pdfeatures[0]==chromosome_label]
		number_overlapping_elements = len(pd_chr.index)
# 		previous_end_points = 0
# 		previous_start_points = 0
		for index,point in pd_chr.iterrows():
			current_height_point=point[2]-point[1]
			current_height_point_scaled=((y_end-y_start)/chr_size)*current_height_point
			point_scaled_start_coord=y_start+((y_end-y_start)/chr_size)*point[1]
			point_scaled_end_coord=y_start+((y_end-y_start)/chr_size)*point[2]
			point_color=colorDictAnnotation[point[column]]
			r=Rectangle(((x_end+annotation_coord),point_scaled_start_coord),.01,current_height_point_scaled,color=point_color)
			ax.add_patch(r)
# 			if ((point_scaled_start_coord.iloc[0]<previous_end_points)|(point_scaled_end_coord.iloc[0]>previous_start_points)):#-0.0002,+0.0002
# 				annotation_coord=annotation_coord+(DIM*.015)
# 			else:
# 				annotation_coord=(DIM*.015)
# 			previous_end_points=point_scaled_end_coord.iloc[0]
# 			previous_start_points=point_scaled_start_coord.iloc[0]
			annotation_coord =annotation_coord+(DIM*.015) # dense
	ax.text(x_start,y_end-(DIM*0.07),'{0} - {1} elements'.format(chromosome_label,number_overlapping_elements))
	collectPatches=[]
	for key,value in colorDictAnnotation.iteritems():
		patch=patches.Patch(color=value,label=key)
		collectPatches.append(patch)
	plt.legend(handles=collectPatches)
	ax.set_yticklabels([])
	ax.set_xticklabels([])

def plot_karyogram(annotationFiles,regions,column,annotationcolumn,genome,spread):
	# get the chromosomes and their sizes to plot, ordered
	sortedChr = select_chromosomes(genome)

	# plot params
	sns.set_style('white')
	fig, ax = plt.subplots()
	
	# set color palette to as many unique column values
	uniqueValsKaryogram = len(regions[annotationcolumn].unique())
	colorPaletteKaryogram = sns.color_palette(palette="Greys",n_colors=uniqueValsKaryogram,desat=.9)
	colorDictKaryogram = dict(zip(regions[annotationcolumn].unique(),colorPaletteKaryogram))
	colorDictAnnotation = unique_annotation_element_dictionary(annotationFiles,column)

	# pdf params
	pp = PdfPages('Karyogram.pdf')
	if spread:
		unique_chromosomes = regions[0].unique()
		for chromosome in unique_chromosomes:
			length = sortedChr.loc[sortedChr['chromosome'] == chromosome]
			chromosome_regions = regions.loc[regions[0] == chromosome]
			set_plot_params_spread(annotationFiles,column,annotationcolumn,chromosome_regions,length,chromosome,ax,colorDictKaryogram,colorDictAnnotation)
			sns.despine(left=True,bottom=True)
			plt.savefig(pp, format='pdf')
			ax.clear()
	else:
		# get the first chromosome for each page, will use to scale others to relative size
		length_chr_one = sortedChr.loc[sortedChr['chromosome'] == 'chr1']
		
		# iterate through the chromosomes, and the patches to plot on the chromosomes
		set_plot_params_dense(annotationFiles,regions,column,annotationcolumn,sortedChr[:12],length_chr_one,ax,colorDictKaryogram,colorDictAnnotation)
		sns.despine(left=True,bottom=True)
		plt.savefig(pp, format='pdf')
		ax.clear()
	
		# iterate through the chromosomes, and the patches to plot on the chromosomes
		set_plot_params_dense(annotationFiles,regions,column,annotationcolumn,sortedChr[12:],length_chr_one,ax,colorDictKaryogram,colorDictAnnotation)
		sns.despine(left=True,bottom=True)
		pp.savefig()
	pp.close()

def main():
	# get args
	args = get_args()

	pdRegions = get_data(args.karyogramdata)
	pdGenome = get_data(args.genome)
	annotationFiles = [line.strip() for line in args.file]
	
	plot_karyogram(annotationFiles,pdRegions,args.colorcolumn,args.annotationcolorcolumn,pdGenome,args.spread)

if __name__ == "__main__":
	main()
