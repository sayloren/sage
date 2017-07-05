"""
Script to perform methylation analyses

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
import tempfile
from ElementLibrary import eachFileProcess
from ElementLibrary import getFeatures
from ElementLibrary import bedtoolToPanda
from ElementLibrary import saveBedTool
from ElementLibrary import pandaToBedtool
from ElementLibrary import reverseComplement

# get methylation positions for all methylation files
def methPositions(mFiles,rangeFeatures,num,uce,inuce,methThresh):
	outMeth = []
	outTable = []
	for methName in mFiles:
		methFeatures = eachFileProcess(methName)
		pdmethThresh = methThreshold(methFeatures,methThresh)
		methPosition = methIntersect(rangeFeatures,pdmethThresh,num,uce,inuce)
		methCpGPos, methCpGNeg, methTable = methMatch(methPosition,rangeFeatures,num) # process the table...
		methCpGPos.rename(columns=lambda x: x.replace('methFreq','{0}_CpGMethylationPos'.format(methName)),inplace=True)
		methCpGNeg.rename(columns=lambda x: x.replace('methFreq','{0}_CpGMethylationNeg'.format(methName)),inplace=True)
		pdmethFreq = methylationFreq(methPosition,num)
		pdmethFreq.columns = ['{0}_Percentage'.format(methName),'{0}_Frequency'.format(methName),'{0}_Coverage'.format(methName)]
		methTable.columns = ['{0}_PosCContext'.format(methName),'{0}_PosMethContext'.format(methName),'{0}_NegCContext'.format(methName),'{0}_NegMethContext'.format(methName)]
		frames1 = [pdmethFreq,methCpGPos,methCpGNeg]
		pdMethEx = pd.concat(frames1,axis=1)
		outMeth.append(pdMethEx)
		outTable.append(methTable)
	pdMeth = pd.concat(outMeth,axis=1)
	pdTable = pd.concat(outTable,axis=1)
	return pdMeth, pdTable

# Threshold the uncapped coverage
def methThreshold(methFeatures,methThresh):
	pdmethFeatures = bedtoolToPanda(methFeatures)
	pdmethThresh = (pdmethFeatures[pdmethFeatures.loc[:,3] >= methThresh])
	btmethThresh = pandaToBedtool(pdmethThresh)
	return btmethThresh

# Intersect the methylation file with the elements; expects 0 based meth, on cytosine
def methIntersect(rangeFeatures,methFeature,num,uce,inuce):
	#savePanda(rangeFeatures[['sBoundary','sEdge','sCenter','eCenter','eEdge','eBoundary']],"Range.txt")
	methSBoundary = methFeature.intersect(rangeFeatures[['chr','sBoundary','start','id']].values.tolist(),wb=True,wa=True)
	if len(methSBoundary) != 0:
		pdmethSBoundary = bedtoolToPanda(methSBoundary)
		pdmethSBoundary['int'] = ((num - uce)/2)
		pdmethSBoundary.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sBoundary','start','id','int']
		pdmethSBoundary['methLoc'] = pdmethSBoundary['int'].astype(int)-(pdmethSBoundary['mstart'].astype(int)-pdmethSBoundary['sBoundary'].astype(int)) # or just count?
		outpdmethSBoundary = pdmethSBoundary[['id','methPer','methLoc','methCov']]
	else:
		outpdmethSBoundary = None

	methSEdge = methFeature.intersect(rangeFeatures[['chr','start','sEdge','id']].values.tolist(),wb=True,wa=True)
	if len(methSEdge) != 0:
		pdmethSEdge = bedtoolToPanda(methSEdge)
		pdmethSEdge['int'] = (((num - uce)/2) + inuce)
		pdmethSEdge.columns = ['mchr','mstart','mstop','methCov','methPer','chr','start','sEdge','id','int']
		pdmethSEdge['methLoc'] = pdmethSEdge['int'].astype(int)-(pdmethSEdge['mstart'].astype(int)-pdmethSEdge['start'].astype(int))
		outpdmethSEdge = pdmethSEdge[['id','methPer','methLoc','methCov']]
	else:
		outpdmethSEdge = None

	methMiddle = methFeature.intersect(rangeFeatures[['chr','sCenter','eCenter','id']].values.tolist(),wb=True,wa=True)
	if len(methMiddle) != 0:
		pdmethFeature = bedtoolToPanda(methMiddle)
		pdmethFeature['int'] = (((num - uce)/2) + (uce - inuce))
		pdmethFeature.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sCenter','eCenter','id','int']
		pdmethFeature['methLoc'] = pdmethFeature['int'].astype(int)-(pdmethFeature['mstart'].astype(int)-pdmethFeature['sCenter'].astype(int))
		outpdmethFeature = pdmethFeature[['id','methPer','methLoc','methCov']]
	else:
		outpdmethFeature = None

	methEEdge = methFeature.intersect(rangeFeatures[['chr','eEdge','end','id']].values.tolist(),wb=True,wa=True)
	if len(methEEdge) != 0:
		pdmethEEdge = bedtoolToPanda(methEEdge)
		pdmethEEdge['int'] = (((num - uce)/2) + uce)
		pdmethEEdge.columns = ['mchr','mstart','mstop','methCov','methPer','chr','eEdge','end','id','int']
		pdmethEEdge['methLoc'] = pdmethEEdge['int'].astype(int)-(pdmethEEdge['mstart'].astype(int)-pdmethEEdge['eEdge'].astype(int))
		outpdmethEEdge = pdmethEEdge[['id','methPer','methLoc','methCov']]
	else:
		outpdmethEEdge = None

	methEBoundary = methFeature.intersect(rangeFeatures[['chr','end','eBoundary','id']].values.tolist(),wb=True,wa=True)
	if len(methEBoundary) != 0:
		pdmethEBoundary = bedtoolToPanda(methEBoundary)
		pdmethEBoundary['int'] = num
		pdmethEBoundary.columns = ['mchr','mstart','mstop','methCov','methPer','chr','end','eBoundary','id','int']#,'methDir'
		pdmethEBoundary['methLoc'] = pdmethEBoundary['int'].astype(int)-(pdmethEBoundary['mstart'].astype(int)-pdmethEBoundary['end'].astype(int))# mstop-stop+startint?
		outpdmethEBoundary = pdmethEBoundary[['id','methPer','methLoc','methCov']]
	else:
		outpdmethEBoundary = None

	methList = [outpdmethSBoundary,outpdmethSEdge,outpdmethFeature,outpdmethEEdge,outpdmethEBoundary]#,
	concatMeth = pd.concat(methList)
	sortMeth = concatMeth.sort_values(['methLoc'],ascending=True)
	return sortMeth

# match cpg and methylation locations
def methMatch(sortMeth,rangeFeatures,num):
	stringDF = rangeFeatures[['id','combineString']]
	groupMeth = methID(sortMeth)
	groupMeth.set_index(keys='id',inplace=True,drop=True)
	groupPosCpG = cpgWindow(stringDF,num,'C')
	groupNegCpG = cpgWindow(stringDF,num,'G')
	stringDF.set_index(keys='id',inplace=True,drop=True)
	methCPos, methTablePos = collapseMethCpG(groupMeth,groupPosCpG,stringDF,num)
	methTablePos.columns = ['PosCContext','PosMethContext']
	methCNeg, methTableNeg = collapseMethCpG(groupMeth,groupNegCpG,stringDF,num)
	methTableNeg.columns = ['NegCContext','NegMethContext']
	methTableNeg.index = methTableNeg.index.map(lambda x: reverseComplement(x))
	methTableNeg.index = methTableNeg.index.str[::-1]
	frames = [methTablePos,methTableNeg]
	methTable = pd.concat(frames,axis=1)
	methTable.fillna('0',inplace=True)
	return methCPos, methCNeg, methTable

# collect methylation locations by id
def methID(sortMeth):
	# but may need to incorporate the other data, (coverage and %) in some way?
	groupMeth = sortMeth.groupby(['id'])['methLoc'].apply(list).to_frame(name = 'groupMeth').reset_index()
	return groupMeth

# return locations of cpgs attached to the id they came from
def cpgWindow(stringDF,num,stFeature):
	out = []
	for element, id in zip(stringDF['combineString'],stringDF['id']):
		index = 0
		while index < len(element): # or num
			index = element.find(stFeature,index)# or 'GC'
			if index == -1:
				break
			outtemp = []
			outtemp.append(id)
			outtemp.append(index)
			index += len(stFeature)
			out.append(outtemp)
	pdout = pd.DataFrame(out)
	pdout.columns = ['id','cpg']
	groupCpG = pdout.groupby(['id'])['cpg'].apply(list).to_frame(name = 'groupCpG').reset_index()
	groupCpG.set_index(keys='id',inplace=True,drop=True)
	return groupCpG

# collapse the methylation and cpg frames, make into data frame by position, make tables to cpg contexts
def collapseMethCpG(groupMeth,groupCpG,stringDF,num):
	frames = [groupMeth,groupCpG,stringDF]
	groupCat = pd.concat(frames,axis=1)
	groupCat.fillna('',inplace=True)#[-1]
	# Use the locations of the c's and the locations of methylation to return exactly where methylation is happening and in what context
	groupCat['groupOverlap'] = groupCat.apply(lambda row:[i for i in row['groupMeth'] if i in row['groupCpG']],axis=1) # instances of methylation and Cs
	groupCat['contextCpG'] = groupCat.apply(lambda row: [row['combineString'][i:i+2] for i in row['groupCpG']],axis=1) # the cpgs at very start,i-1:i+2
	groupCat['contextOverlap'] = groupCat.apply(lambda row: [row['combineString'][i:i+2] for i in row['groupOverlap']],axis=1) # the cpgs at very start,i-1:i+2
# 	groupCat['methNoOv'] = groupCat.apply(lambda row: [sorted(set(row['groupMeth']) - set(row['groupCpG']))],axis=1)
	groupCat['ContextDict'] = groupCat.apply(lambda row: [dict(zip(row['groupOverlap'],row['contextOverlap']))],axis=1)
	Context = groupCat['ContextDict'].apply(pd.Series).stack().tolist()
	filterContext = filter(None, Context)
	splitContext = {}
	for d in filterContext:
		for k, v in d.iteritems():
			splitContext.setdefault(v,[]).append(k)
	outMethC = []
	outColNames = []
	for k, v in splitContext.items():
		colName = 'methFreq{0}'.format(k)
		outColNames.append(colName)
		methC = methylatedCpGFreq(v,num)
		outMethC.append(methC)
	pdSepContext = pd.concat(outMethC,axis=1)
	pdSepContext.columns = outColNames
	Methylated = groupCat['groupOverlap'].apply(pd.Series).stack().tolist()
	methC = methylatedCpGFreq(Methylated,num)
	methCat = pd.concat([methC,pdSepContext],axis=1)
	methCat = methCat.rename(columns = {'methFreq':'methFreqTotal'})
	ContextCpG = groupCat['contextCpG'].apply(pd.Series).stack().tolist()
	ContextCpG = filter(lambda s: len(s) > 1,ContextCpG) # cpgs at the end will be filtered out, 2
	contextCCount = pd.Series(ContextCpG).value_counts()
	ContextOverlap = groupCat['contextOverlap'].apply(pd.Series).stack().tolist()
	ContextOverlap = filter(lambda s: len(s) > 1, ContextOverlap) # cpgs at the end will be filtered out,2
	contextOverlapCount = pd.Series(ContextOverlap).value_counts()
	frames = [contextCCount,contextOverlapCount]
	contextTable = pd.concat(frames,axis=1)
	contextTable.columns = ['CContext','MethContext']
	return methCat, contextTable

# collect methylated (cpg) frequencies by position of elements
def methylatedCpGFreq(Methylated,num): # concatMeth from methIntersect
	if not Methylated:
		methCpG = pd.DataFrame(0,index=np.arange(num),columns=['methFreq'])
	else:
		new_index = range(0,num)
		pdMeth = pd.DataFrame(Methylated).astype(int)
		pdMeth.columns = ['methLoc']
		pdMeth['methFreq'] = pdMeth.groupby(['methLoc'])['methLoc'].transform('count')
		dupMeth = pdMeth.sort_values(['methLoc'],ascending=True).drop_duplicates(['methLoc','methFreq']) # returns highest value of methylation
		methCpG = dupMeth.set_index('methLoc').reindex(new_index,fill_value=0)
		methCpG.index.name = None
	return methCpG

# collect methylation frequencies by position of elements
def methylationFreq(methPosition,num): # concatMeth from methIntersect
	new_index = range(0,num)
	subsetMeth = methPosition[['methLoc','methPer','methCov']]
	subsetMeth['methFreq'] = subsetMeth.groupby(['methLoc'])['methPer'].transform('count')
	MethPer = subsetMeth[['methLoc','methPer','methFreq']]
	MethCov = subsetMeth[['methLoc','methCov','methFreq']]
	sortPer = MethPer.sort_values(['methLoc','methPer','methFreq'],ascending=True)
	sortCov = MethCov.sort_values(['methLoc','methCov','methFreq'],ascending=True)
	dupPer = sortPer.drop_duplicates(['methLoc','methFreq'],keep='last')# returns highest value
	dupCov = sortCov.drop_duplicates(['methLoc','methFreq'],keep='last')# returns highest value
	indexPer = dupPer.set_index('methLoc').reindex(new_index,fill_value=0)
	indexCov = dupCov.set_index('methLoc').reindex(new_index,fill_value=0)
	indexPer.index.name = None
	indexCov.index.name = None
	intPer = indexPer.astype(int)
	intCov = indexCov.astype(int)
	frames = [intPer,intCov]
	intMeth = pd.concat(frames,axis=1)
	outMeth = intMeth.loc[:,~intMeth.columns.duplicated()]
	return outMeth

def main(mFiles,rangeFeatures,num,uce,inuce,methThresh):
	pdMeth, pdTable = methPositions(mFiles,rangeFeatures,num,uce,inuce,methThresh)
	return pdMeth, pdTable

if __name__ == "__main__":
	main()
