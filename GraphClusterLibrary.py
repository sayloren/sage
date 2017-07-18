"""
Script to graph the clustermap analysis, and return clustered info
Otherwise the settings mess with the other graphs

Wren Saylor
July 12 2017

To Do:
tissue x uce
tissue x meth/c
uce x meth/c
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
from GraphMethLibrary import methIndex

# Make some graphs for fangs
def graphCluster(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine):

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
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Cluster_{0}.pdf'.format(fileName))

	# Use the row_colors to color those with similar SD?
# 	lut = dict(zip(species.unique(), "rbg"))
# 	row_colors = species.map(lut)
	heatmap0 = sns.clustermap(ATgroup,cmap='RdPu',vmin=0,vmax=100,xticklabels=100,col_cluster=False)#,row_colors=row_colors
	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap0.ax_heatmap.set_yticks([]))
	plt.setp(heatmap0.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap0.ax_heatmap.set_ylabel('{0} UCEs'.format(len(ATgroup.index)),size=8))
	plt.setp(heatmap0.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap0.ax_heatmap.set_title('Mean AT Content per UCE',size=12))
	
	sns.despine()
	pp.savefig()

	# Various combinations to plot on heatmaps
	# Frequency x (tissue, id)
	FreqPlusTis,FreqMinusTis = methIndex(pdMeth,'tissue','methFreq',num)
	FreqPlusID,FreqMinusID = methIndex(pdMeth,'id','methFreq',num)
	
	# Make heatmap for # methylation on pos strand (Frequency)
	ylabels1 = FreqPlusTis.index
	heatmap1 = sns.clustermap(FreqPlusTis,cmap='RdPu',xticklabels=100,col_cluster=False,yticklabels=ylabels1)#cbar_ax=cbar5_ax,vmin=0,vmax=5
	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap1.ax_heatmap.set_ylabel('Sample',size=8))
	plt.setp(heatmap1.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap1.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap1.ax_heatmap.set_yticks(np.arange(FreqPlusTis.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap1.ax_heatmap.set_yticklabels(ylabels1,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap1.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))
	
	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	ylabels2 = FreqMinusTis.index
	heatmap2 = sns.clustermap(FreqMinusTis,cmap='RdPu',xticklabels=100,col_cluster=False,yticklabels=ylabels2)#cbar_ax=cbar5_ax#,vmin=0,vmax=5
	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap2.ax_heatmap.set_ylabel('Sample',size=8))
	plt.setp(heatmap2.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap2.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap2.ax_heatmap.set_yticks(np.arange(FreqMinusTis.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap2.ax_heatmap.set_yticklabels(ylabels2,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap2.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap3 = sns.clustermap(FreqPlusID,cmap='RdPu',xticklabels=100,col_cluster=False)#,vmin=0,vmax=5
	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap3.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqPlusID.index)),size=8))
	plt.setp(heatmap3.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap3.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap3.ax_heatmap.set_yticks([]))
	plt.setp(heatmap3.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap4 = sns.clustermap(FreqMinusID,cmap='RdPu',xticklabels=100,col_cluster=False)#,vmin=0,vmax=5
	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap4.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqMinusID.index)),size=8))
	plt.setp(heatmap4.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap4.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap4.ax_heatmap.set_yticks([]))
	plt.setp(heatmap4.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))
# 	#http://seaborn.pydata.org/examples/timeseries_bootstrapped.html

	sns.despine()
	pp.savefig()
	pp.close()

def main(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine):
	graphCluster(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine)

if __name__ == "__main__":
	main()