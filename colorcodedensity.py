"""
Script to graph a density score (pairing) with color coded bars on the underlyingraph

Jumana Alhaj Abed
December 2017

Copyright 2017 Harvard University, Wu Lab

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

import pybedtools as pbt
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.gridspec as gridspec

# set command line arguments
def get_args():
	# file lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("-d","--densityfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?
	parser.add_argument("-c","--chromatincolorsfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?
	return parser.parse_args()

# set some global variables
def set_global_variables(args):
	global density
	density = agrs.densityfile
	global chromatin
	chromatin = args.chromatincolorsfile

def main():
	args = get_args()
	set_global_variables(args)
	
	# 1) read in the bed files
		# I defiantly have a couple variations of this around
	
	# 2) change to pandas
		# this should be around too
	
	# 3) graph the density file
		# something like plt.plot()
		# threshold line - plt.axhline(y,linewidth=.05,linestyle='dashed')
	
	# 4) change the colors on the graph background
		# https://stackoverflow.com/questions/9957637/how-can-i-set-the-background-color-on-specific-areas-of-a-pyplot-figure
		# plt.axvspan(i, i+.5, facecolor='b', alpha=0.5)
	
	# 5) do 2) and 3), but threshold before plotting
		# probably will help if you add na's where you don't want to graph
		# you use gridspec to manipulate where the subplots
		# https://matplotlib.org/users/gridspec.html
	
	# 6) save pdf
		# I like using pdfpages
		# https://matplotlib.org/xkcd/examples/pylab_examples/multipage_pdf.html
	
if __name__ == "__main__":
	main()