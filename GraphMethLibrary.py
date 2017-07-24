"""
Script to graph the methylation analyses

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
import seaborn as sns

# Transform the Frequency, Percentage and Coverage data into graphable data frames
def methIndex(dataframe,yItem,zItem,num): 
	# x item is methLoc, y item is either tissue or id, z item is coverage, percentage, or frequency
	new_index = range(0,num)
	
	# Separate by strand
	PlusMeth = dataframe.loc[dataframe['Cytosine'] == 'C']
	MinusMeth = dataframe.loc[dataframe['Cytosine'] == 'G']
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[['methLoc',yItem,zItem]]
	MinussubMeth = MinusMeth[['methLoc',yItem,zItem]]

	# Sort ascending, in order to only use the highest value with keep = last
	PlussortMeth = PlussubMeth.sort_values(['methLoc'],ascending=True).drop_duplicates(['methLoc',yItem,zItem],keep='last')
	MinussortMeth = MinussubMeth.sort_values(['methLoc'],ascending=True).drop_duplicates(['methLoc',yItem,zItem],keep='last')

	# Pivot the data frame so that each tissue/cell type is a column
	PluspivotMeth = pd.pivot_table(PlussortMeth,index='methLoc',columns=yItem,values=zItem)
	MinuspivotMeth = pd.pivot_table(MinussortMeth,index='methLoc',columns=yItem,values=zItem)
	
	PluspivotMeth.columns.name = None
	MinuspivotMeth.columns.name = None
	
	# Give new index, using the methLocations
	PlusindexMeth = PluspivotMeth.reindex(new_index,fill_value=0)
	MinusindexMeth = MinuspivotMeth.reindex(new_index,fill_value=0)

	# Fill in missing index with 0s
	PlusindexMeth.fillna('0',inplace=True)
	MinusindexMeth.fillna('0',inplace=True)
	
	# Remove the index column name
	PlusindexMeth.index.name = None
	MinusindexMeth.index.name = None
	
	# Transpose the data frame for easy input into the heatamp
	PlustransMeth = PlusindexMeth.T
	MinustransMeth = MinusindexMeth.T
	
	PlustransMeth = PlustransMeth[PlustransMeth.columns].astype(float)
	MinustransMeth = MinustransMeth[MinustransMeth.columns].astype(float)
	
	print 'Converted {0} by {1} into data frame'.format(yItem,zItem)
	
	return PlustransMeth, MinustransMeth

# Make Methylation graphs
def graphMeth(slidingWinDF,pdMeth,fileName,num,uce,inuce,window):
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)
	sns.set_style('ticks')
	info = str(fileName) + ', '+ str(len(slidingWinDF)) + ' - ' "UCES"
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Methylation_{0}.pdf'.format(fileName))

	# Various combinations to plot on heatmaps
	# x Tissue
	FreqPlusTis,FreqMinusTis = methIndex(pdMeth,'tissue','methFreq',num)
	PerPlusTis,PerMinusTis = methIndex(pdMeth,'tissue','methPer',num)
	CovPlusTis,CovMinusTis = methIndex(pdMeth,'tissue','methCov',num)

	# x ID
	FreqPlusID,FreqMinusID = methIndex(pdMeth,'id','methFreq',num)
	PerPlusID,PerMinusID = methIndex(pdMeth,'id','methPer',num)
	CovPlusID,CovMinusID = methIndex(pdMeth,'id','methCov',num)
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	
	# Make heatmap for methylation coverage on pos strand, capped at 1000
	# Should make the bottom threshold a different color
	ax0 = plt.subplot(gs[0,:])
	heatmap0 = sns.heatmap(CovPlusTis,cmap='RdPu',ax=ax0,vmin=0,vmax=100,xticklabels=100)
	ax0.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax0.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax0.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax0.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax0.set_ylabel('Sample',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.tick_params(labelsize=8)
	ylabels0 = CovPlusTis.index
	ax0.set_yticklabels(ylabels0,minor=False,rotation=0)
	ax0.set_yticks(np.arange(CovPlusTis.shape[0]) + 0.5, minor=False)
	ax0.set_title('Methylation Coverage for Plus Strand',size=8)

	# Make heatmap for methylation coverage on neg strand, capped at 1000
	# Should make the bottom threshold a different color
	ax1 = plt.subplot(gs[1,:])
	heatmap1 = sns.heatmap(CovMinusTis,cmap='RdPu',ax=ax1,vmin=0,vmax=100,xticklabels=100)
	ax1.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax1.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax1.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax1.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax1.set_ylabel('Sample',size=8)
	ax1.set_xlabel('Position',size=6)
	ax1.tick_params(labelsize=8)
	ylabels1 = CovMinusTis.index
	ax1.set_yticklabels(ylabels1,minor=False,rotation=0)
	ax1.set_yticks(np.arange(CovMinusTis.shape[0]) + 0.5, minor=False)
	ax1.set_title('Methylation Coverage for Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Make heatmap for # methylation on pos strand (Frequency)
	ax2 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap2 = sns.heatmap(FreqPlusTis,cmap='RdPu',ax=ax2,xticklabels=100)#cbar_ax=cbar5_ax,vmin=0,vmax=5
	ax2.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax2.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax2.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax2.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax2.set_ylabel('Sample',size=8)
	ax2.set_xlabel('Position',size=6)
	ax2.tick_params(labelsize=8)
	ylabels2 = FreqPlusTis.index
	ax2.set_yticklabels(ylabels2,minor=False,rotation=0)
	ax2.set_yticks(np.arange(FreqPlusTis.shape[0]) + 0.5, minor=False)
	ax2.set_title('Methylation Frequency on Plus Strand',size=8)

	# Make heatmap for # methylation on pos strand (Frequency)
	ax3 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap3 = sns.heatmap(FreqMinusTis,cmap='RdPu',ax=ax3,xticklabels=100)#cbar_ax=cbar5_ax,vmin=0,vmax=5
	ax3.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax3.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax3.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax3.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax3.set_ylabel('Sample',size=8)
	ax3.set_xlabel('Position',size=6)
	ax3.tick_params(labelsize=8)
	ylabels3 = FreqMinusTis.index
	ax3.set_yticklabels(ylabels3,minor=False,rotation=0)
	ax3.set_yticks(np.arange(FreqMinusTis.shape[0]) + 0.5, minor=False)
	ax3.set_title('Methylation Frequency on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Make heatmap for % methylation on pos strand
	# Should make the bottom threshold a different color
	ax4 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap4 = sns.heatmap(PerPlusTis,cmap='RdPu',ax=ax4,vmin=0,vmax=100,xticklabels=100)
	ax4.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax4.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax4.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax4.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax4.set_ylabel('Sample',size=8)
	ax4.set_xlabel('Position',size=6)
	ax4.tick_params(labelsize=8)
	ylabels4 = PerPlusTis.index
	ax4.set_yticklabels(ylabels4,minor=False,rotation=0)
	ax4.set_yticks(np.arange(PerPlusTis.shape[0]) + 0.5, minor=False)
	ax4.set_title('Methylation Percentage on Plus Strand',size=8)

	# Make heatmap for % methylation on neg strand
	# Should make the bottom threshold a different color
	ax5 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap5 = sns.heatmap(PerMinusTis,cmap='RdPu',ax=ax5,vmin=0,vmax=100,xticklabels=100)
	ax5.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax5.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax5.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax5.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax5.set_ylabel('Sample',size=8)
	ax5.set_xlabel('Position',size=6)
	ax5.tick_params(labelsize=8)
	ylabels5 = PerMinusTis.index
	ax5.set_yticklabels(ylabels5,minor=False,rotation=0)
	ax5.set_yticks(np.arange(PerMinusTis.shape[0]) + 0.5, minor=False)
	ax5.set_title('Methylation Percentage on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Make heatmap for methylation coverage on pos strand, capped at 1000
	# Should make the bottom threshold a different color
	ax6 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap6 = sns.heatmap(CovPlusID,cmap='RdPu',ax=ax6,xticklabels=100,vmin=0,vmax=100)
	ax6.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax6.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax6.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax6.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax6.set_ylabel('{0} UCEs'.format(len(CovPlusID.index)),size=8)
	ax6.set_xlabel('Position',size=6)
	ax6.tick_params(labelsize=8)
	ax6.set_yticklabels([])
	ax6.set_title('Methylation Coverage on Plus Strand',size=8)

	# Make heatmap for methylation coverage on neg strand, capped at 1000
	# Should make the bottom threshold a different color
	ax7 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap7 = sns.heatmap(CovMinusID,cmap='RdPu',ax=ax7,xticklabels=100,vmin=0,vmax=100)
	ax7.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax7.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax7.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax7.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax7.set_ylabel('{0} UCEs'.format(len(CovMinusID.index)),size=8)
	ax7.set_xlabel('Position',size=6)
	ax7.tick_params(labelsize=8)
	ax7.set_yticklabels([])
	ax7.set_title('Methylation Coverage on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Make heatmap for # methylation on pos strand (Frequency)
	ax8 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap8 = sns.heatmap(FreqPlusID,cmap='RdPu',ax=ax8,xticklabels=100)#,vmin=0,vmax=5
	ax8.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax8.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax8.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax8.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax8.set_ylabel('{0} UCEs'.format(len(FreqPlusID.index)),size=8)
	ax8.set_xlabel('Position',size=6)
	ax8.tick_params(labelsize=8)
	ax8.set_yticklabels([])
# 	ax8.set_yticks(np.arange(FreqPlusID.shape[0]) + 0.5, minor=False)
	ax8.set_title('Methylation Frequency on Plus Strand',size=8)

	# Make heatmap for # methylation on neg strand (Frequency)
	ax9 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap9 = sns.heatmap(FreqMinusID,cmap='RdPu',ax=ax9,xticklabels=100)#,vmin=0,vmax=5
	ax9.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax9.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax9.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax9.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax9.set_ylabel('{0} UCEs'.format(len(FreqMinusID.index)),size=8)
	ax9.set_xlabel('Position',size=6)
	ax9.tick_params(labelsize=8)
	ax9.set_yticklabels([])
	ax9.set_title('Methylation Frequency on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Make heatmap for % methylation on pos strand
	# Should make the bottom threshold a different color
	ax10 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap10 = sns.heatmap(PerPlusID,cmap='RdPu',ax=ax10,vmin=0,vmax=100,xticklabels=100)
	ax10.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax10.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax10.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax10.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax10.set_ylabel('{0} UCEs'.format(len(PerPlusID.index)),size=8)
	ax10.set_xlabel('Position',size=6)
	ax10.tick_params(labelsize=8)
	ax10.set_yticklabels([])
	ax10.set_title('Methylatoin Percentage on Plus Strand',size=8)

	# Make heatmap for % methylation on neg strand
	# Should make the bottom threshold a different color
	ax11 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap11 = sns.heatmap(PerMinusID,cmap='RdPu',ax=ax11,vmin=0,vmax=100,xticklabels=100)
	ax11.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax11.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax11.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax11.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax11.set_ylabel('{0} UCEs'.format(len(PerMinusID.index)),size=8)
	ax11.set_xlabel('Position',size=6)
	ax11.tick_params(labelsize=8)
	ax11.set_yticklabels([])
	ax11.set_title('Methylatoin Percentage on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	pp.close()

def main(slidingWinDF,pdMeth,fileName,num,uce,inuce,window):
	graphMeth(slidingWinDF,pdMeth,fileName,num,uce,inuce,window)

if __name__ == "__main__":
	main()