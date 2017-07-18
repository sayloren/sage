"""
Script to graph the Fang analyses

Wren Saylor
July 5 2017
"""

import argparse
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
import scipy.stats as ss
from scipy.stats import mstats
import seaborn as sns

# Make some graphs for fangs
def graphFang(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine):

	# Parameters that all graphs will use
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)

	# Get mean and standard deviation for AT
	ATNames = [names.index(i) for i in names if 'A' in i or 'T' in i]
	ATDataFrames = [slidingWinDF[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()
	ATstd = ATgroup.std()
	
	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1])
	gs.update(hspace=.8) # setting the space between the graphs
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
	
	# Plot the mean AT content with a std of 1
	StartMean = ATgroup.loc[:,(((num-uce)/2)-halfwindow-window):(((num-uce)/2)+(inuce-halfwindow))].mean() # the boundary and 50inward
	StopMean = ATgroup.loc[:,(((num-uce)/2)+(uce-inuce-halfwindow)):(((num-uce)/2)+uce-(halfwindow-window))].mean() # the boundary and 50inward
	wilcoxPSRMean = ss.wilcoxon(StartMean,StopMean)
	ax0 = plt.subplot(gs[0])
	ax0.plot(fillX,ATmean,linewidth=1, color='#3e1638',label='AT')
	ax0.fill_between(fillX,ATmean+ATstd,ATmean-ATstd,facecolor = '#63245a',label='',alpha=0.3)
	ax0.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.hlines(y=86,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
	ax0.text(32,85,'11bp sliding window',size=6)
	ax0.text(20,90,'Wilcox Signed Rank P-value {:0.1e}'.format(wilcoxPSRMean[1]),size=6,clip_on=False)
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax0.set_title('Mean AT Content With Standard Deviation',size=8)
	ax0.set_yticks(ax0.get_yticks()[::2])
	plt.xlim(0,num)

	# Plot the std = 1
	ax1 = plt.subplot(gs[1],sharex=ax0)
	ax1.plot(fillX,ATstd,linewidth=1, color='#3e1638',label='AT')
	ax1.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvspan((((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow),facecolor = '#ae3e9e',label='',alpha=0.2)
	ax1.set_yticks(ax1.get_yticks()[::2])
	ax1.set_xlabel('Position',size=6)
	ax1.set_ylabel('SD',size=8)
	ax1.set_title('Standard Deviation',size=8)
	plt.setp(ax1.get_xticklabels(), visible=True)
	ax1.legend(loc=0,fontsize=5,labelspacing=0.05)
	
	# Significances tests for SD populations
	uceRegion = ATgroup.loc[:,(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)].std()
	bothStream = ATgroup.iloc[:,np.r_[0:(((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow):(num-window)]].std()
	# Kruskal-Wallis test
# 	kruskalSD = mstats.kruskalwallis(bothStream,uceRegion)
# 	ax2.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
	ax2 = plt.subplot(gs[2])
	ax2.hist(uceRegion,35,linewidth=0.3, color='#ae3e9e',label='UCE',alpha=0.5)
	ax2.hist(bothStream,35,linewidth=0.3, color='#ae3e66',label='Surrounding Regions',alpha=0.5)
	ax2.set_yticks(ax2.get_yticks()[::2])
	ax2.set_ylabel('Frequency',size=8)
	ax2.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax2.set_xlabel('Standard Deviation Value',size=8)

	sns.despine()
	plt.savefig(pp, format='pdf')
	
	# Separate out those with only a single nucleotide search
	SingleNames = [names.index(i) for i in names if len(i) == 1]
	SingleNamesVal = [names[i] for i in SingleNames]
	SingleDataFrames = [slidingWinDF[i] for i in SingleNames]
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	sns.set_palette("husl")
	
	# Nucleotides
	ax3 = plt.subplot(gs[0,:],sharex=ax0)
	for dfNuc,lNuc in zip(SingleDataFrames,SingleNamesVal):
		ax3.plot(fillX,dfNuc.mean(),linewidth=1,label=lNuc)
	ax3.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.set_yticks(ax3.get_yticks()[::2])
	ax3.set_ylabel('% Nucleotide Content',size=8)
	ax3.set_xlabel('Position',size=6)
	ax3.set_title('Mean Nucleotide Content',size=8)
	ax3.legend(loc=0,fontsize=5,labelspacing=0.05)
	
	# Plot SD for Nucleotides
	# Kruskal-Wallis test
# 	kruskalSD = mstats.kruskalwallis(uceRegionCpA,uceRegionCpT,uceRegionCpC,uceRegionCpG)
# 	ax2.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
# 	https://stackoverflow.com/questions/35089422/two-seaborn-distplots-one-same-axis
	ax13 = plt.subplot(gs[1,:])
	for dfNuc,lNuc in zip(SingleDataFrames,SingleNamesVal):
		elRegion = dfNuc.loc[:,(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)].std()
# 		elFlank = dfNuc.iloc[:,np.r_[0:(((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow):(num-window)]].std()
		ax13.hist(elRegion,35,linewidth=0.3,alpha=0.5,label='{0}-{1}'.format(lNuc,'element'))
# 		ax13.hist(elFlank,35,linewidth=0.3,alpha=0.5,label='{0}-{1}'.format(lNuc,'flank'))
	ax13.set_yticks(ax13.get_yticks()[::2])
	ax13.set_ylabel('Frequency',size=8)
	ax13.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax13.set_xlabel('Standard Deviation Value',size=8)
	ax13.set_title('Standard Deviation for Nucleotides within Element',size=8)

	# For any search stings that are dinucleotides
	if any(len(i) == 2 for i in names):
		
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)
		sns.set_palette("husl")
		
		# Separate out those with only a double nucleotide search
		DoubleNames = [names.index(i) for i in names if len(i) == 2]
		DoubleNamesVal = [names[i] for i in DoubleNames]
		DoubleDataFrames = [slidingWinDF[i] for i in DoubleNames]
		
		# Plot the CpN
		# Might still want to return the actual CpN location for how many are methylated
		ax14 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(DoubleDataFrames,DoubleNamesVal):
			ax14.plot(fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax14.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
		ax14.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
		ax14.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
		ax14.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
		ax14.set_title('Mean Nucleotide String Content',size=8)
		ax14.set_ylabel('% Nucleotide String Content',size=8)
		ax14.set_xlabel('Position',size=6)
		ax14.set_yticks(ax14.get_yticks()[::2])
		ax14.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for CpNs
# 		Kruskal-Wallis test
# 		kruskalSD = mstats.kruskalwallis(uceRegionCpA,uceRegionCpT,uceRegionCpC,uceRegionCpG)
# 		ax15.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
		ax15 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(DoubleDataFrames,DoubleNamesVal):
			elRegion = dfNuc.loc[:,(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)].std()
			# elFlank = dfNuc.iloc[:,np.r_[0:(((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow):(num-window)]].std()
			ax15.hist(elRegion,25,linewidth=0.3,label=lNuc,alpha=0.5)
		ax15.set_yticks(ax15.get_yticks()[::2])
		ax15.set_ylabel('Frequency',size=8)
		ax15.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax15.set_xlabel('Standard Deviation Value',size=8)
		ax15.set_title('Standard Deviation for Nucleotide String within Element',size=8)

	# For any search strings 3 bp long
	if any(len(i) == 3 for i in names):
	
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)
		sns.set_palette("husl")
		
		# Separate out those with only a multiple nucleotide search
		TriNames = [names.index(i) for i in names if len(i) == 3]
		TriNamesVal = [names[i] for i in TriNames]
		TriDataFrames = [slidingWinDF[i] for i in TriNames]

		# Plot the MultiNucleotide Sequences
		# Might still want to return the actual CpN location for how many are methylated
		ax16 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(TriDataFrames,TriNamesVal):
			ax16.plot(fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax16.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
		ax16.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
		ax16.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
		ax16.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
		ax16.set_title('Mean Trinucleotide Content',size=8)
		ax16.set_ylabel('% Trinucleotide Content',size=8)
		ax16.set_xlabel('Position',size=6)
		ax16.set_yticks(ax16.get_yticks()[::2])
		ax16.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for MultiNucleotide Sequences
# 		Kruskal-Wallis test
# 		kruskalSD = mstats.kruskalwallis(uceRegionCpA,uceRegionCpT,uceRegionCpC,uceRegionCpG)
# 		ax15.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
		ax17 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(TriDataFrames,TriNamesVal):
			elRegion = dfNuc.loc[:,(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)].std()
			#elFlank = dfNuc.iloc[:,np.r_[0:(((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow):(num-window)]].std()
			ax17.hist(elRegion,20,linewidth=0.3,label=lNuc,alpha=0.5)
		ax17.set_yticks(ax17.get_yticks()[::2])
		ax17.set_ylabel('Frequency',size=8)
		ax17.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax17.set_xlabel('Standard Deviation Value',size=8)
		ax17.set_title('Standard Deviation for Trinucleotides within Element',size=8)

	# For any search strings longer than 3
	if any(len(i) > 3 for i in names):
	
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)
		sns.set_palette("husl")
		
		# Separate out those with only a multiple nucleotide search
		MultiNames = [names.index(i) for i in names if len(i) > 3]
		MultiNamesVal = [names[i] for i in MultiNames]
		MultiDataFrames = [slidingWinDF[i] for i in MultiNames]

		# Plot the MultiNucleotide Sequences
		# Might still want to return the actual CpN location for how many are methylated
		ax18 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(MultiDataFrames,MultiNamesVal):
			ax18.plot(fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax18.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
		ax18.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
		ax18.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
		ax18.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
		ax18.set_title('Mean Polynucleotide Content',size=8)
		ax18.set_ylabel('% Polynucleotide Content',size=8)
		ax18.set_xlabel('Position',size=6)
		ax18.set_yticks(ax18.get_yticks()[::2])
		ax18.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for MultiNucleotide Sequences
# 		Kruskal-Wallis test
# 		kruskalSD = mstats.kruskalwallis(uceRegionCpA,uceRegionCpT,uceRegionCpC,uceRegionCpG)
# 		ax15.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
		ax19 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(MultiDataFrames,MultiNamesVal):
			elRegion = dfNuc.loc[:,(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)].std()
			#elFlank = dfNuc.iloc[:,np.r_[0:(((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow):(num-window)]].std()
			ax19.hist(elRegion,20,linewidth=0.3,label=lNuc,alpha=0.5)
		ax19.set_yticks(ax19.get_yticks()[::2])
		ax19.set_ylabel('Frequency',size=8)
		ax19.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax19.set_xlabel('Standard Deviation Value',size=8)
		ax19.set_title('Standard Deviation for Polynucleotides within Element',size=8)

	sns.despine()
	pp.savefig()
	pp.close()

# Collect each UCEs second derivative
def behaviorUCE(fillX,pdWindow):
	secondderUCE = []
	for index, row in pdWindow.iterrows():
		f = splrep(fillX,row,k=5,s=11)
		smoothMean = splev(fillX,f)
		secondDer = splev(fillX,f,der=2)
		secondDer[0:window] = 0 # small edge effect
		secondDer[-window:] = 0 # small edge effect
		secondderUCE.append(secondDer)
	pdSecderUCE = pd.DataFrame(secondderUCE)
	return pdSecderUCE

def main(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine):
	graphFang(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine)

if __name__ == "__main__":
	main()