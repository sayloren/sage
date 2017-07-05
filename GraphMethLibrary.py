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
from numpy import sin, linspace, pi
import numdifftools as nd
import numpy as np
import numpy.ma as ma
import pandas as pd
import re
import scipy as sp
import scipy.fftpack
from scipy.interpolate import splrep, splev
from scipy import signal as ssignal
import scipy.stats as ss
from scipy.stats import mstats
import seaborn as sns

# make some graphs!
def graphMeth(pdWindow,pdMeth,pdTable,fileName,num,uce,inuce,window):
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)
	sns.set_style('ticks')
	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1]) # for the ratios of the graphs against each other
	gs.update(hspace=.8) # setting the space between the graphs
	info = str(fileName) + ', '+ str(len(pdWindow)) + ' - ' "UCES"
	plt.suptitle(info,fontsize=10)
	mean = pdWindow.mean()
	#maX = np.ma.masked_where(240<fillX<250 & 340<fillX<350, fillX) # tried to make a mast for gapped data
	pp = PdfPages('Methylation_{0}.pdf'.format(fileName))

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Table Methylation Contexts
	pdTable1 = (pdTable[pdTable.columns[pdTable.columns.str.contains('PosCContext',case=False)]])
	pdTable2 = (pdTable[pdTable.columns[pdTable.columns.str.contains('PosMethContext',case=False)]])
	pdTable1 = pdTable1.mean(axis=1)
	outTable1 = pd.concat([pdTable1,pdTable2],axis=1)
	outTable1.rename(columns={0:"Count"},inplace=True)
	pdTable3 = (pdTable[pdTable.columns[pdTable.columns.str.contains('NegCContext',case=False)]])
	pdTable4 = (pdTable[pdTable.columns[pdTable.columns.str.contains('NegMethContext',case=False)]])
	pdTable3 = pdTable3.mean(axis=1)
	outTable2 = pd.concat([pdTable3,pdTable4],axis=1)
	outTable2.rename(columns={0:"Count"},inplace=True)
# 	pdTable1['PerMeth'] = (pdTable1['CContextSum']/pdTable1['MethContextSum'])* 100

	ax18 = plt.subplot(gs[0,:])
	ax18.set_frame_on(False)
	ax18.set_yticks([])
	ax18.set_xticks([])
	ylabels3 = outTable1.columns.str.replace('.bed_PosMethContext','')
	MethTable1 = ax18.table(cellText=outTable1.values,rowLabels=outTable1.index,colLabels=ylabels3,cellLoc='center',rowLoc='center',loc='center')
	ax18.set_title('Plus Strand Methylation',size=8)
	MethTable1.set_fontsize(8)

	ax19 = plt.subplot(gs[1,:])
	ax19.set_frame_on(False)
	ax19.set_yticks([])
	ax19.set_xticks([])
	ylabels4 = outTable2.columns.str.replace('.bed_NegMethContext','')
	MethTable2 = ax19.table(cellText=outTable2.values,rowLabels=outTable2.index,colLabels=ylabels4,cellLoc='center',rowLoc='center',loc='center')
	ax19.set_title('Minus Strand Methylation',size=8)
	MethTable2.set_fontsize(8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	
	# Those that are methylated
	pdCpGPosTotal = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationPosTotal',case=False)]])
	pdCpGNegTotal = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationNegTotal',case=False)]])
	pdCpGPosCpG = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationPosCG',case=False)]])
	pdCpGNegGpC = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationNegGC',case=False)]])
	pdCpGPosCpC = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationPosCC',case=False)]])
	pdCpGNegGpG = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationNegGG',case=False)]])
	pdCpGPosCpA = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationPosCA',case=False)]])
	pdCpGNegGpT = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationNegGT',case=False)]])
	pdCpGPosCpT = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationPosCT',case=False)]])
	pdCpGNegGpA = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('CpGMethylationNegGA',case=False)]])
	
	# Methylation Data Intersections
	pdMethPer = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('Percentage',case=False)]])
	pdMethNum = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('Frequency',case=False)]])
	pdMethCov = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('Coverage',case=False)]])
	# Transposes
	TPer = pdMethPer.T
	TNum = pdMethNum.T
	TCov = pdMethCov.T
	TPos = pdCpGPosTotal.T
	TNeg = pdCpGNegTotal.T
	TPosCG = pdCpGPosCpG.T
	TNegGC = pdCpGNegGpC.T
	TPosCC = pdCpGPosCpC.T
	TNegGG = pdCpGNegGpG.T
	TPosCA = pdCpGPosCpA.T
	TNegGT = pdCpGNegGpT.T
	TPosCT = pdCpGPosCpT.T
	TNegGA = pdCpGNegGpA.T

	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	gs.update(hspace=.5)

		# Make heatmap for % methylation (Percentage)
	sns.set_style('ticks')
	ax26 = plt.subplot(gs[0,:])
	ylabels11 = TPer.index.str.replace('.bed_Percentage','')
	heatmap11 = sns.clustermap(TPer,cmap='RdPu',vmin=0,vmax=100,xticklabels=100,yticklabels=ylabels11,col_cluster=False)#ax=ax26,metric="correlation"
	#http://seaborn.pydata.org/examples/timeseries_bootstrapped.html
	ax26.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax26.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax26.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax26.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	plt.setp(heatmap11.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap11.ax_heatmap.set_yticks(np.arange(TPer.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap11.ax_heatmap.set_yticklabels(ylabels11,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap11.ax_heatmap.set_ylabel('Sample',size=8))
	plt.setp(heatmap11.ax_heatmap.set_xlabel('Position',size=6))
	plt.setp(heatmap11.ax_heatmap.tick_params(labelsize=8))
# 	ax26.set_yticklabels(minor=False,rotation=0)#ylabels11,
# 	ax26.set_yticks(np.arange(TPer.shape[0]) + 0.5, minor=False)
	ax26.set_title('Methylation Percentage',size=8)

	sns.despine()
	pp.savefig()
	plt.clf()
	plt.cla()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Make heatmap for methylation coverage, capped at 1000
	ax27 = plt.subplot(gs[0,:],sharex=ax26)
	heatmap12 = sns.heatmap(TCov,cmap='RdPu',ax=ax27,vmin=0,vmax=100,xticklabels=100)
	ax27.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax27.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax27.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax27.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax27.set_ylabel('Sample',size=8)
	ax27.set_xlabel('Position',size=6)
	ax27.tick_params(labelsize=8)
	ylabels12 = TCov.index.str.replace('.bed_Coverage','')
	ax27.set_yticklabels(ylabels12,minor=False,rotation=0)
	ax27.set_yticks(np.arange(TCov.shape[0]) + 0.5, minor=False)
	ax27.set_title('Methylation Coverage',size=8)
	
	# Heat map for interaction of % and #
	
	# Make heatmap for # methylation (Frequency)
	ax28 = plt.subplot(gs[1,:],sharex=ax26)
	heatmap13 = sns.heatmap(TNum,cmap='RdPu',ax=ax28,vmin=0,vmax=5,xticklabels=100)#cbar_ax=cbar5_ax
	ax28.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax28.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax28.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax28.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax28.set_ylabel('Sample',size=8)
	ax28.set_xlabel('Position',size=6)
	ax28.tick_params(labelsize=8)
	ylabels13 = TNum.index.str.replace('.bed_Frequency','')
	ax28.set_yticklabels(ylabels13,minor=False,rotation=0)#
	ax28.set_yticks(np.arange(TNum.shape[0]) + 0.5, minor=False)
	ax28.set_title('Methylation Frequency',size=8)
	#https://stackoverflow.com/questions/35869150/creating-a-diverging-color-palette-with-a-midrange-instead-of-a-midpoint
	#http://seaborn.pydata.org/examples/structured_heatmap.html
	
	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	
	#['#ffffff','#f5eef1','#af7994','#a56887','#9b587a','#8b4f6d','#7c4661','#6c3d55','#5d3449','#4d2c3d','#3e2330','#2e1a24','#1f1118','#0f080c','#000000']
	# Make heatmap for Positive Strand Total Cytosine Methylation
	ax16 = plt.subplot(gs[0,:],sharex=ax26)
	heatmap1 = sns.heatmap(TPos,cmap=ListedColormap(['#f5eef1','#9b587a','#7c4661','#5d3449','#3e2330','#1f1118']),vmin=0,vmax=5,ax=ax16,xticklabels=100)#'RdPu'
	ax16.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax16.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax16.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax16.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax16.set_ylabel('Sample',size=8)
	ax16.set_xlabel('Position',size=6)
	ax16.tick_params(labelsize=8)
	ylabels1 = TPos.index.str.replace('.bed_CpGMethylationPosTotal','')
	ax16.set_yticklabels(ylabels1,minor=False,rotation=0)
	ax16.set_yticks(np.arange(TPos.shape[0]) + 0.5, minor=False)
	ax16.set_title('Total Cytosine Methylation on Plus Strand',size=8)

	# Make heatmap for Negitive Strand Total Cytosine Methylation
	ax17 = plt.subplot(gs[1,:],sharex=ax26)
	heatmap2 = sns.heatmap(TNeg,cmap=ListedColormap(['#f5eef1','#9b587a','#7c4661','#5d3449','#3e2330','#1f1118']),vmin=0,vmax=5,ax=ax17,xticklabels=100)#'RdPu'
	ax17.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax17.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax17.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax17.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax17.set_ylabel('Sample',size=8)
	ax17.set_xlabel('Position',size=6)
	ax17.tick_params(labelsize=8)
	ylabels2 = TNeg.index.str.replace('.bed_CpGMethylationNegTotal','')
	ax17.set_yticklabels(ylabels2,minor=False,rotation=0)
	ax17.set_yticks(np.arange(TNeg.shape[0]) + 0.5, minor=False)
	ax17.set_title('Total Cytosine Methylation on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# '#fbebe8','#c13918','#962c12','#6b200d','#401308','#150602'
	# Make heatmap for Positive Strand CpG Methylation
	ax18 = plt.subplot(gs[0,:],sharex=ax26)
	heatmap3 = sns.heatmap(TPosCG,cmap=ListedColormap(['#fbebe8','#c13918','#962c12','#6b200d','#401308','#150602']),ax=ax18,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap,'RdPu'
	ax18.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax18.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax18.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax18.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax18.set_ylabel('Sample',size=8)
	ax18.set_xlabel('Position',size=6)
	ax18.tick_params(labelsize=8)
	ylabels3 = TPosCG.index.str.replace('.bed_CpGMethylationPosCG','')
	ax18.set_yticklabels(ylabels3,minor=False,rotation=0)
	ax18.set_yticks(np.arange(TPosCG.shape[0]) + 0.5, minor=False)
	ax18.set_title('CpG Methylation on Plus Strand',size=8)

	# Make heatmap for Negitive Strand CpG Methylation
	ax19 = plt.subplot(gs[1,:],sharex=ax26)
	heatmap4 = sns.heatmap(TNegGC,cmap=ListedColormap(['#fbebe8','#c13918','#962c12','#6b200d','#401308','#150602']),ax=ax19,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax19.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax19.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax19.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax19.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax19.set_ylabel('Sample',size=8)
	ax19.set_xlabel('Position',size=6)
	ax19.tick_params(labelsize=8)
	ylabels4 = TNegGC.index.str.replace('.bed_CpGMethylationNegGC','')
	ax19.set_yticklabels(ylabels4,minor=False,rotation=0)
	ax19.set_yticks(np.arange(TNegGC.shape[0]) + 0.5, minor=False)
	ax19.set_title('CpG Methylation on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# '#fbe8f7','#c118a0','#96127c','#6b0d59','#400835','#150211'
	# Make heatmap for Positive Strand CpC Methylation
	ax20 = plt.subplot(gs[0,:],sharex=ax26)
	heatmap5 = sns.heatmap(TPosCC,cmap=ListedColormap(['#fbe8f7','#c118a0','#96127c','#6b0d59','#400835','#150211']),ax=ax20,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax20.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax20.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax20.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax20.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax20.set_ylabel('Sample',size=8)
	ax20.set_xlabel('Position',size=6)
	ax20.tick_params(labelsize=8)
	ylabels5 = TPosCC.index.str.replace('.bed_CpGMethylationPosCC','')
	ax20.set_yticklabels(ylabels5,minor=False,rotation=0)
	ax20.set_yticks(np.arange(TPosCC.shape[0]) + 0.5, minor=False)
	ax20.set_title('CpC Methylation on Plus Strand',size=8)

	# Make heatmap for Negitive Strand CpC Methylation
	ax21 = plt.subplot(gs[1,:],sharex=ax26)
	heatmap6 = sns.heatmap(TNegGG,cmap=ListedColormap(['#fbe8f7','#c118a0','#96127c','#6b0d59','#400835','#150211']),ax=ax21,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax21.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax21.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax21.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax21.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax21.set_ylabel('Sample',size=8)
	ax21.set_xlabel('Position',size=6)
	ax21.tick_params(labelsize=8)
	ylabels6 = TNegGG.index.str.replace('.bed_CpGMethylationNegGG','')
	ax21.set_yticklabels(ylabels6,minor=False,rotation=0)
	ax21.set_yticks(np.arange(TNegGG.shape[0]) + 0.5, minor=False)
	ax21.set_title('CpC Methylation on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# '#f5e8fb','#8e18c1','#6e1296','#4f0d6b','#2f0840','#0f0215'
	# Make heatmap for Positive Strand CpA Methylation
	ax22 = plt.subplot(gs[0,:],sharex=ax26)
	heatmap7 = sns.heatmap(TPosCA,cmap=ListedColormap(['#f5e8fb','#8e18c1','#6e1296','#4f0d6b','#2f0840','#0f0215']),ax=ax22,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax22.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax22.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax22.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax22.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax22.set_ylabel('Sample',size=8)
	ax22.set_xlabel('Position',size=6)
	ax22.tick_params(labelsize=8)
	ylabels7 = TPosCA.index.str.replace('.bed_CpGMethylationPosCA','')
	ax22.set_yticklabels(ylabels7,minor=False,rotation=0)
	ax22.set_yticks(np.arange(TPosCA.shape[0]) + 0.5, minor=False)
	ax22.set_title('CpA Methylation on Plus Strand',size=8)

	# Make heatmap for Negitive Strand CpA Methylation
	ax23 = plt.subplot(gs[1,:],sharex=ax26)
	heatmap8 = sns.heatmap(TNegGT,cmap=ListedColormap(['#f5e8fb','#8e18c1','#6e1296','#4f0d6b','#2f0840','#0f0215']),ax=ax23,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax23.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax23.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax23.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax23.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax23.set_ylabel('Sample',size=8)
	ax23.set_xlabel('Position',size=6)
	ax23.tick_params(labelsize=8)
	ylabels8 = TNegGT.index.str.replace('.bed_CpGMethylationNegGT','')
	ax23.set_yticklabels(ylabels8,minor=False,rotation=0)
	ax23.set_yticks(np.arange(TNegGT.shape[0]) + 0.5, minor=False)
	ax23.set_title('CpA Methylation on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# '#fbe8ed','#c1184b','#96123a','#6b0d2a','#400819','#150208'
	# Make heatmap for Positive Strand CpT Methylation
	ax24 = plt.subplot(gs[0,:],sharex=ax26)
	heatmap9 = sns.heatmap(TPosCT,cmap=ListedColormap(['#fbe8ed','#c1184b','#96123a','#6b0d2a','#400819','#150208']),ax=ax24,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax24.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax24.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax24.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax24.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax24.set_ylabel('Sample',size=8)
	ax24.set_xlabel('Position',size=6)
	ax24.tick_params(labelsize=8)
	ylabels9 = TPosCT.index.str.replace('.bed_CpGMethylationPosCT','')
	ax24.set_yticklabels(ylabels9,minor=False,rotation=0)
	ax24.set_yticks(np.arange(TPosCT.shape[0]) + 0.5, minor=False)
	ax24.set_title('CpT Methylation on Plus Strand',size=8)

	# Make heatmap for Negitive Strand CpT Methylation
	ax25 = plt.subplot(gs[1,:],sharex=ax26)
	heatmap10 = sns.heatmap(TNegGA,cmap=ListedColormap(['#fbe8ed','#c1184b','#96123a','#6b0d2a','#400819','#150208']),ax=ax25,vmin=0,vmax=5,xticklabels=100)#col_cluster=False sns.clustermap
	ax25.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax25.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax25.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax25.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax25.set_ylabel('Sample',size=8)
	ax25.set_xlabel('Position',size=6)
	ax25.tick_params(labelsize=8)
	ylabels10 = TNegGA.index.str.replace('.bed_CpGMethylationNegGA','')
	ax25.set_yticklabels(ylabels10,minor=False,rotation=0)
	ax25.set_yticks(np.arange(TNegGA.shape[0]) + 0.5, minor=False)
	ax25.set_title('CpT Methylation on Minus Strand',size=8)

	sns.despine()
	pp.savefig()
	plt.clf()
	plt.cla()
	pp.close()

def main(pdWindow,pdMeth,pdTable,fileName,num,uce,inuce,window):
	graphMeth(pdWindow,pdMeth,pdTable,fileName,num,uce,inuce,window)

if __name__ == "__main__":
	main()