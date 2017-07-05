"""
Script to perform RC sorting

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
from FangsLibrary import slidingWindow
from ElementLibrary import eachFileProcess
from MethylationLibrary import methPositions
from MethylationLibrary import methThreshold
from MethylationLibrary import methIntersect
from MethylationLibrary import methMatch
from MethylationLibrary import methylationFreq

def dirLine(directionFeatures,mFiles,num,uce,inuce,window,methThresh):
	negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
	outComp, outCpA, outCpT, outCpG, outCpC, compA, compT, compG, compC, compMo, outnegMeth, outnegTable  = [], [], [], [], [], [], [], [], [], [], [],[]
	for element in negStr['reverseComplement']:
		negFeature, negCpAfeature,negCpTfeature,negCpGfeature,negCpCfeature,negAfeature,negTfeature,negGfeature,negCfeature,negMofeature = slidingWindow(element,num,uce,inuce,window)
		outComp.append(negFeature)
		outCpA.append(negCpAfeature)
		outCpT.append(negCpTfeature)
		outCpG.append(negCpGfeature)
		outCpC.append(negCpCfeature)
		compA.append(negAfeature)
		compT.append(negTfeature)
		compG.append(negGfeature)
		compC.append(negCfeature)
		compMo.append(negMofeature)
	for methName in mFiles:
		methFeatures = eachFileProcess(methName)
		pdmethThresh = methThreshold(methFeatures,methThresh)
		methPosition = methIntersect(negStr,pdmethThresh,num,uce,inuce)
		methCpGPos, methCpGNeg, methTable = methMatch(methPosition,directionFeatures,num)
		methCpGPos.rename(columns=lambda x: x.replace('methFreq','{0}_CpGMethylationPos'.format(methName)),inplace=True)
		methCpGNeg.rename(columns=lambda x: x.replace('methFreq','{0}_CpGMethylationNeg'.format(methName)),inplace=True)
		methFreq = methylationFreq(methPosition,num)
		methFreq.columns = ['{0}_Percentage'.format(methName),'{0}_Frequency'.format(methName),'{0}_Coverage'.format(methName)]
		methTable.columns = ['{0}_PosCContext'.format(methName),'{0}_PosMethContext'.format(methName),'{0}_NegCContext'.format(methName),'{0}_NegMethContext'.format(methName)]
		frames = [methFreq,methCpGPos,methCpGNeg]
		pdMethEx = pd.concat(frames,axis=1)
		outnegMeth.append(pdMethEx)
		outnegTable.append(methTable)
	posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
	for element in posStr['combineString']:
		posFeature, posCpAfeature,posCpTfeature,posCpGfeature,posCpCfeature,posAfeature, posTfeature, posGfeature, posCfeature,posMofeature = slidingWindow(element,num,uce,inuce,window)
		outComp.append(posFeature)
		outCpA.append(posCpAfeature)
		outCpT.append(posCpTfeature)
		outCpG.append(posCpGfeature)
		outCpC.append(posCpCfeature)
		compA.append(posAfeature)
		compT.append(posTfeature)
		compG.append(posGfeature)
		compC.append(posCfeature)
		compMo.append(posMofeature)
	compWindow,compWinCpA,compWinCpT,compWinCpG,compWinCpC,compWinA,compWinT,compWinG,compWinC,compWinMo = pd.DataFrame(outComp),pd.DataFrame(outCpA),pd.DataFrame(outCpT),pd.DataFrame(outCpG),pd.DataFrame(outCpC),pd.DataFrame(compA),pd.DataFrame(compT),pd.DataFrame(compG),pd.DataFrame(compC),pd.DataFrame(compMo)
	pdposMeth, pdposTable = methPositions(mFiles,posStr,num,uce,inuce,methThresh)
	pdnegMeth = pd.concat(outnegMeth,axis=1)
	pdnegTable = pd.concat(outnegTable,axis=1)
	pdnegInMeth = pdnegMeth.iloc[::-1] # reverse order 
	pdnegInsetMeth = pdnegInMeth.reset_index(drop=True) # make new index, don't make old one into column
	pdcompMeth = pd.concat([pdnegInsetMeth,pdposMeth],axis=1)
	pdcompTable = pd.concat([pdnegTable,pdposTable],axis=1)
	pdcompTable = pdcompTable.astype(int)
	pdgroupTable = pdcompTable.groupby(pdcompTable.columns,axis=1).sum()
	pdcompMethPer = pdcompMeth.loc[:,pdcompMeth.columns.str.contains('Percentage',case=False)]
	pdgroupMethPer = pdcompMethPer.groupby(pdcompMethPer.columns, axis=1).max()
	pdcompMethNum = pdcompMeth.loc[:,pdcompMeth.columns.str.contains('Frequency',case=False)]
	pdgroupMethNum = pdcompMethNum.groupby(pdcompMethNum.columns, axis=1).sum()
	pdcompMethCov = pdcompMeth.loc[:,pdcompMeth.columns.str.contains('Coverage',case=False)]
	pdgroupMethCov = pdcompMethCov.groupby(pdcompMethCov.columns, axis=1).max()
	pdcompMethPos = pdcompMeth.loc[:,pdcompMeth.columns.str.contains('CpGMethylationPos',case=False)]
	pdgroupMethPos = pdcompMethPos.groupby(pdcompMethPos.columns, axis=1).sum()
	pdcompMethNeg = pdcompMeth.loc[:,pdcompMeth.columns.str.contains('CpGMethylationNeg',case=False)]
	pdgroupMethNeg = pdcompMethNeg.groupby(pdcompMethNeg.columns, axis=1).sum()
	pdgroupMeth = pd.concat([pdgroupMethPer,pdgroupMethNum,pdgroupMethCov,pdgroupMethPos,pdgroupMethNeg],axis=1)
	return pdgroupMeth,pdgroupTable,compWindow,compWinCpA,compWinCpT,compWinCpG,compWinCpC,compWinA,compWinT,compWinG,compWinC,compWinMo

def main(directionFeatures,binDir,mFiles,num,uce,inuce,window,methThresh):
	pdgroupMeth,pdgroupTable,compWindow,compWinCpA,compWinCpT,compWinCpG,compWinCpC,compWinA,compWinT,compWinG,compWinC,compWinMo = dirLine(directionFeatures,mFiles,num,uce,inuce,window,methThresh)
	return pdgroupMeth,pdgroupTable,compWindow,compWinCpA,compWinCpT,compWinCpG,compWinCpC,compWinA,compWinT,compWinG,compWinC,compWinMo

if __name__ == "__main__":
	main()