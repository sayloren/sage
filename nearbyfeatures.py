"""

Script to look at how close features are between two bed files, or a UCE file and a list of bedfiles

Wren Saylor
Jan 2018

Files to run:
-UCEs
-Domains (not pooled)
-Genes

To Do:
-Stats per bed file features (number, size, genome coverage)
-Stats per overlaps (min/max, average, std) 
-Where do overlaps occur

"""
import argparse
import pandas as pd
import pybedtools as pbt

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=str,help='the primary element file')
	parser.add_argument("-s","--secondaryfeatures",type=argparse.FileType('rU'),help="a list of files with the secondary features to query")
	parser.add_argument("-t","--tertiaryfeatures",type=argparse.FileType('rU'),help="a list of files with the tertiary features to query")
	return parser.parse_args()




def main():




if __name__ == "__main__":
	main()