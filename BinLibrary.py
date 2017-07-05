"""

Script to determine how the probability of getting boundaries with the same
value changes as the bin size of the element boundaries changes

Wren Saylor
June 20 2017

"""

import argparse
import itertools
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb
import seaborn as sns

def graphComb(n,k):
	#(n+range(k)-1) / (i)
	# order not important and repetition allowed
	krange = np.arange(k)
	combinations = comb(n,krange, repetition=True)
	sns.set_style('ticks')
	plt.plot(krange,combinations)
	sns.despine()
	plt.savefig('Combinations.pdf')

def main(n,k):
	graphComb(n,k)

if __name__ == "__main__":
	main()