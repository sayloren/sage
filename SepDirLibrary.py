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

# separate by directionality
def dirLine(directionFeatures,mFiles,num,uce,inuce,window,methThresh):
	dirWindow, dirWinCpG, dirWinA, dirWinT, dirWinG, dirWinC, dirWinMo = dataframeWindow(dirStr,num,uce,inuce,window)
	pddirMeth, pddirTable = methPositions(mFiles,dirStr,num,uce,inuce,methThresh)
	return pddirMeth, pddirTable, dirWindow, dirWinCpG, dirWinA, dirWinT, dirWinG, dirWinC, dirWinMo

def main(directionFeatures,mFiles,num,uce,inuce,window,methThresh):
	pddirMeth,pddirTable,dirWindow,dirWinCpG,dirWinA,dirWinT,dirWinG,dirWinC,dirWinMo = dirLine(directionFeatures,mFiles,num,uce,inuce,window,methThresh)
	return pddirMeth,pddirTable,dirWindow,dirWinCpG,dirWinA,dirWinT,dirWinG,dirWinC,dirWinMo

if __name__ == "__main__":
	main()