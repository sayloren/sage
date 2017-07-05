""""""

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
import tempfile

# run append to sliding window and return pandas data frames
def dataframeWindow(rangeFeatures,num,uce,inuce,window):
	outWindow,outCpA,outCpT,outCpG,outCpC,outA,outT,outG,outC, outMo = appendWindow(rangeFeatures,num,uce,inuce,window)
	pdWindow,pdCpA,pdCpT,pdCpG,pdCpC,pdA,pdT,pdG,pdC,pdMo = pd.DataFrame(outWindow),pd.DataFrame(outCpA),pd.DataFrame(outCpT),pd.DataFrame(outCpG),pd.DataFrame(outCpC),pd.DataFrame(outA),pd.DataFrame(outT),pd.DataFrame(outG),pd.DataFrame(outC),pd.DataFrame(outMo)
	return pdWindow, pdCpA, pdCpT,pdCpG,pdCpC, pdA, pdT, pdG, pdC, pdMo

# append results from sliding window
def appendWindow(rangeFeatures,num,uce,inuce,window):
	outWindow = []
	outCpA = []
	outCpT = []
	outCpG = []
	outCpC = []
	outA = []
	outT = []
	outG = []
	outC = []
	outMo = []
	for element in rangeFeatures['combineString']:
		outFeature, winCpA, winCpT, winCpG, winCpC, winA, winT, winG, winC, winMo = slidingWindow(element,num,uce,inuce,window)
		outWindow.append(outFeature)
		outCpA.append(winCpA)
		outCpT.append(winCpT)
		outCpG.append(winCpG)
		outCpC.append(winCpC)
		outA.append(winA)
		outT.append(winT)
		outG.append(winG)
		outC.append(winC)
		outMo.append(winMo)
	return outWindow, outCpA, outCpT, outCpG, outCpC, outA, outT, outG, outC, outMo

# sliding window
def slidingWindow(element,num,uce,inuce,window):
	winFeatures = []
	winCpA = []
	winCpT = []
	winCpG = []
	winCpC = []
	winA = []
	winT = []
	winG = []
	winC = []
	winMotif = []
	n = num #600 # len(element) # take out hard coding
	start, end = 0, window
	while end < n:
		current = element[start:end]
		#print num, float(len(element)), float(len(current)), start, end, current
		percentageAT = eval('100 * float(current.count("A") + current.count("T"))/ float(len(current))')
		percentageCpA = eval('100 * float(current.count("CA")) / float(len(current))')
		percentageCpT = eval('100 * float(current.count("CT")) / float(len(current))')
		percentageCpG = eval('100 * float(current.count("CG")) / float(len(current))')
		percentageCpC = eval('100 * float(current.count("CC")) / float(len(current))')
		perA = eval('100 * float(current.count("A"))/ float(len(current))')
		perT = eval('100 * float(current.count("T"))/ float(len(current))')
		perG = eval('100 * float(current.count("G"))/ float(len(current))')
		perC = eval('100 * float(current.count("C"))/ float(len(current))')
		perMo = eval('100 * float(current.count("ATTAAT")) / float(len(current))')
		winFeatures.append(percentageAT)
		winCpA.append(percentageCpA)
		winCpT.append(percentageCpT)
		winCpG.append(percentageCpG)
		winCpC.append(percentageCpC)
		winA.append(perA)
		winT.append(perT)
		winG.append(perG)
		winC.append(perC)
		winMotif.append(perMo)
		start, end = start + 1, end + 1
	return winFeatures, winCpA, winCpT, winCpG, winCpC, winA, winT, winG, winC, winMotif

def main(rangeFeatures,num,uce,inuce,window):
	pdWindow, pdCpA, pdCpT, pdCpG, pdCpC, pdA, pdT, pdG, pdC, pdMo = dataframeWindow(rangeFeatures,num,uce,inuce,window)
	return pdWindow, pdCpA, pdCpT, pdCpG, pdCpC, pdA, pdT, pdG, pdC, pdMo 
	
if __name__ == "__main__":
	main()