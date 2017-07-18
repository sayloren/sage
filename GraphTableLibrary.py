"""
Script to return tables 

Wren Sayor
July 13 2017

To Do:
Need to calculate the total cpn across the entire genome!!!
Put each tissue type on separate page
"""

import argparse
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Make Methylation graphs
def graphTable(slidingWinDF,pdTable,fileName,num,uce,inuce,window):
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)
	sns.set_style('ticks')
	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1])
	gs.update(hspace=.8)
	info = str(fileName) + ', '+ str(len(slidingWinDF)) + ' - ' "UCES"
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Table_{0}.pdf'.format(fileName))

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	# Table Methylation Contexts
# 	pdTable1 = (pdTable[pdTable.columns[pdTable.columns.str.contains('PosCContext',case=False)]])
# 	pdTable2 = (pdTable[pdTable.columns[pdTable.columns.str.contains('PosMethContext',case=False)]])
# 	pdTable1 = pdTable1.mean(axis=1)
# 	outTable1 = pd.concat([pdTable1,pdTable2],axis=1)
# 	outTable1.rename(columns={0:"Count"},inplace=True)
# 	pdTable3 = (pdTable[pdTable.columns[pdTable.columns.str.contains('NegCContext',case=False)]])
# 	pdTable4 = (pdTable[pdTable.columns[pdTable.columns.str.contains('NegMethContext',case=False)]])
# 	pdTable3 = pdTable3.mean(axis=1)
# 	outTable2 = pd.concat([pdTable3,pdTable4],axis=1)
# 	outTable2.rename(columns={0:"Count"},inplace=True)
# 	pdTable1['PerMeth'] = (pdTable1['CContextSum']/pdTable1['MethContextSum'])* 100

	ax0 = plt.subplot(gs[0,:])
	ax0.set_frame_on(False)
	ax0.set_yticks([])
	ax0.set_xticks([])
# 	ylabels1 = outTable1.columns.str.replace('.bed_PosMethContext','')
# 	MethTable1 = ax0.table(cellText=outTable1.values,rowLabels=outTable1.index,colLabels=ylabels1,cellLoc='center',rowLoc='center',loc='center')
	ax0.set_title('Plus Strand Methylation',size=8)
	MethTable1.set_fontsize(8)

	ax1 = plt.subplot(gs[1,:])
	ax1.set_frame_on(False)
	ax1.set_yticks([])
	ax1.set_xticks([])
# 	ylabels2 = outTable2.columns.str.replace('.bed_NegMethContext','')
# 	MethTable2 = ax1.table(cellText=outTable2.values,rowLabels=outTable2.index,colLabels=ylabels2,cellLoc='center',rowLoc='center',loc='center')
	ax1.set_title('Minus Strand Methylation',size=8)
	MethTable2.set_fontsize(8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	pp.close()

def main(slidingWinDF,pdTable,fileName,num,uce,inuce,window):
	graphTable(slidingWinDF,pdTable,fileName,num,uce,inuce,window)

if __name__ == "__main__":
	main()