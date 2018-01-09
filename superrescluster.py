"""
Script to make some graphs from the vutara clustering data

Wren Saylor
October 2017
"""

import argparse
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pltq
import seaborn as sns
import matplotlib as mpl
import matplotlib.gridspec as gridspec

def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("-f","--file",type=str,help="file from vutara clustering")
# 	parser.add_argument("-c","--clusters",type=int,nargs='*',help="file from vutara clustering")
	return parser.parse_args()

def read_csv_file(file):
	pdfile = pd.read_csv(file)
	return pdfile

def line_graphs(pdcsv):
# 	particlecountdf = pdcsv.loc[:,~pdcsv.columns.str.contains('PCA|centroid')]
	cols = pdcsv.columns.tolist()
	newcols = [i.encode('ascii', 'ignore') for i in cols]
	pdcsv.columns = newcols
	particledf = pdcsv.set_index('cluster',drop=True)
	particledf.index.name=None
	particledf = particledf.rename(columns=lambda x: x.strip())
	domains = [1,2,4,5,7]
	ucedf = particledf[particledf['probe'].isin(domains)]
	
	sns.set_style('ticks')
	sns.set_palette("husl",n_colors=5)
	gs = gridspec.GridSpec(3,1,height_ratios=[1,1,1])
	gs.update(hspace=.8)
	pp = PdfPages('three_graphs.pdf')
# 	plt.figure(figsize=(14,7))
	ax0 = plt.subplot(gs[0,:])
	ax1 = plt.subplot(gs[1,:],sharex=ax0)
	ax2 = plt.subplot(gs[2,:],sharex=ax0)
	
	sns.pointplot(ucedf.index,ucedf['density'],hue=ucedf['probe'],ax=ax0,legend_out=True,size=10,loc=2)
	sns.pointplot(ucedf.index,ucedf['area'],hue=ucedf['probe'],ax=ax1,legend_out=True,size=10,loc=2)
	sns.pointplot(ucedf.index,ucedf['particleCount'],hue=ucedf['probe'],ax=ax2,legend_out=True,size=10,loc=2)
	
	sns.despine()
	pp.savefig()
	pp.close()


def main():
	args = get_args()
	pdcsv = read_csv_file(args.file)
	
	line_graphs(pdcsv)

if __name__ == "__main__":
	main()
