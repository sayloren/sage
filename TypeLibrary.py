"""
Script to separate by type

Wren Saylor
July 5 2017
"""
import argparse
from collections import defaultdict
import itertools
import numpy as np
import pandas as pd
import pybedtools as pbt
from MethylationLibrary import methPositions
from FangsLibrary import dataframeWindow

# do all the analysis for each type
def perType(boolType,fileName,binDir,mFiles,num,uce,inuce,window,methThresh):
	pdWindow, pdTypeCpA,pdTypeCpT,pdTypeCpG,pdTypeCpC, pdTypeA, pdTypeT, pdTypeG, pdTypeC, pdTypeMo = dataframeWindow(boolType,num,uce,inuce,window)
	pdMeth, pdTable = methPositions(mFiles,boolType,num,uce,inuce,methThresh)
	return pdMeth,pdTable,pdWindow,pdTypeCpA,pdTypeCpT,pdTypeCpG,pdTypeCpC,pdTypeA,pdTypeT,pdTypeG,pdTypeC,pdTypeMo

def main(boolType,fileName,binDir,mFiles,num,uce,inuce,window,methThresh):
	pdMeth,pdTable,pdWindow,pdTypeCpA,pdTypeCpT,pdTypeCpG,pdTypeCpC,pdTypeA,pdTypeT,pdTypeG,pdTypeC,pdTypeMo = perType(boolType,fileName,binDir,mFiles,num,uce,inuce,window,methThresh)
	return pdMeth,pdTable,pdWindow,pdTypeCpA,pdTypeCpT,pdTypeCpG,pdTypeCpC,pdTypeA,pdTypeT,pdTypeG,pdTypeC,pdTypeMo

if __name__ == "__main__":
	main()
