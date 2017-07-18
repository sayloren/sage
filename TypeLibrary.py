"""
Script to separate by type

Wren Saylor
July 5 2017
"""
import argparse
import pandas as pd
from MethylationLibrary import compactMeth
from FangsLibrary import compactWindow

# do all the analysis for each type
def perType(boolType,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome):
	typeWindow, typeNames = compactWindow(boolType['combineString'],boolType['id'],num,uce,inuce,window,nucLine)
	pdMeth = compactMeth(mFiles,boolType,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	return pdMeth,typeWindow,typeNames

def main(boolType,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome):
	pdMeth,typeWindow,typeNames = perType(boolType,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
	print 'Completed group sorting for {0} items'.format(len(boolType.index))
	return pdMeth,typeWindow,typeNames

if __name__ == "__main__":
	main()
