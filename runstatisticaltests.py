"""
Script to take in two csv with two columns and run several statistical tests
Returns description of the test/ what data sets it is appropriate for, coefficient and p-value heatmaps

Wren Saylor
January 2018

Guy Requests: QQplot, skeweness
"""

import argparse
import pandas as pd

def get_args():
	parser.add_argument("file", type=str,help='A file with two columns to run statistical tests on')
# 	parser.add_argumnet("-t",'--testoptions',type=str,choices=[],help='statistical test options')
# 	parser.add_argument("-s",'--sectiondata',type=int,nargs='4',help='choose section of two columns to run statistical tests on, in format [col 1 start, col 1 end, col 2 start, col 2 end]')
	return parser.parse_args()

def statistical_test_dataframe():
	kolmogorov_smirnov = ['ks','nonparametric test, testing equality of continuous, one-dimensional probability distributions']
	wilcoxon_signed_rank = ['wsr','nonparamtric test to assess whether population mean ranks differ']
	kruskal_wallis = ['kw','nonparametric method testing whether samples originate from the same distribution, extends mann-whitney u test']
	mann_whitney_wilcoxon = ['mww','nonparametic test for likelyhood that a random value from one sample will be >/< a random value from a second sample']
	students_t = ['stt','parametric test determing signifint difference between data sets']
	spearman_rank_correlation = ['src','nonparametic test for seeing whether the ranks of two variables covary']
	
	
	makedf = pd.Dataframe()
	return makedf

def main():
	args = get_args()
	
	# 1) 
	# print list of stat test and descriptions
	
	# 2)
	# interactive command line to select tests
	
	# 3)
	# run tests
	
	# 4)
	# print tests run with descriptions, heat-maps of coefficients and p values

if __name__ == "__main__":
	main()