"""
Script to separate by directionality

Wren Saylor
July 5 2017
"""
import argparse
import Bio
from Bio import SeqIO
from Bio import Seq
from collections import defaultdict
import itertools
from numpy import sin, linspace, pi
import numdifftools as nd
import numpy as np
import numpy.ma as ma
import pandas as pd
import pybedtools as pbt
import re
import scipy as sp
import scipy.fftpack
from scipy.interpolate import splrep, splev
from scipy import signal as ssignal
import scipy.stats as ss
from scipy.stats import mstats
import tempfile
from MethylationLibrary import methPositions
from FangsLibrary import dataframeWindow

# out put directionality, as inferred by comparing first and last n base pairs from input parameters
def compareN(element,size):
	start = element[:size]
	end = element[-size:]
	perSize = []
	perSize.append(eval('100*int(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*int(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
	# give + - = depending on which side has larger AT content
	if perSize[0] > perSize[1]: outList = '+'
	if perSize[1] > perSize[0]: outList = '-'
	if perSize[1] == perSize[0]: outList = '='
	return outList

# with the results from compareN per each element, evaluate directionality into new column
def evalN(rangeFeatures,fileName,binDir):
	rangeFeatures['compareBoundaries'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],binDir)),axis=1)
	compareEnds = pd.DataFrame(rangeFeatures[['chr','start','end','compareBoundaries']])
	return rangeFeatures

def main(rangeFeatures,fileName,binDir):
	directionFeatures = evalN(rangeFeatures,fileName,binDir)
	return directionFeatures

if __name__ == "__main__":
	main()