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

"""
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "{}"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntplt for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright {yyyy} {name of copyright owner}

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