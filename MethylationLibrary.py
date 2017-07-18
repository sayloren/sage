"""
Script to perform methylation analyses

Wren Saylor
July 5 2017

"""
import argparse
import pandas as pd
from ElementLibrary import eachFileProcess
from ElementLibrary import getFeatures
from ElementLibrary import bedtoolToPanda
from ElementLibrary import saveBedTool
from ElementLibrary import pandaToBedtool
from ElementLibrary import reverseComplement
from ElementLibrary import simpleFasta

# 	Old snippets from previous methylation processing method
# 	df.rename(columns=lambda x: x.replace('string','{0}_string'.format(filename)),inplace=True) # modify column names to be file specific
# 	df.set_index(keys='id',inplace=True,drop=True) # change index to column id
# 	df = df.iloc[::-1] # reverse order of index
# 	df = df.groupby(df.columns,axis=1).sum() # group by column name and sum those with same name
# 	dfstring = df.loc[:,df.columns.str.contains('string',case=False)] # search for column containing string and subset into new df
# 	df['group1'] = df.apply(lambda row:[i for i in row['group2'] if i in row['group3']],axis=1) # get a new column where two columns lists overlap
# 	df['string'] = df.apply(lambda row: [row['longstring'][i:i+2] for i in row['location']],axis=1) # get the local string for a location in the string
# 	df['dictionary'] = df.apply(lambda row: [dict(zip(row['group1'],row['group2']))],axis=1) # make a list a dictionary for another list in column
# 	List = df['list'].apply(pd.Series).stack().tolist() # make a list from a column of lists

def methThreshold(methFeatures,methCovThresh,methPerThresh):
	pdmethFeatures = bedtoolToPanda(methFeatures)
	pdmethCovThresh = (pdmethFeatures[pdmethFeatures.loc[:,3] >= methCovThresh])
	pdmethPerThresh = (pdmethCovThresh[pdmethCovThresh.loc[:,4] >= methPerThresh])
	btmethThresh = pandaToBedtool(pdmethPerThresh)
	print 'Methylation coverage is being thresholded at {0} and percentage at {1}'.format(methCovThresh, methPerThresh)
	return btmethThresh

def methIntersect(rangeFeatures,methFeature,num,uce,inuce,faGenome):
	methSBoundary = methFeature.intersect(rangeFeatures[['chr','sBoundary','sEdge','id']].values.tolist(),wb=True,wa=True)
	if len(methSBoundary) != 0:
		pdmethSBoundary = bedtoolToPanda(methSBoundary)
		pdmethSBoundary['int'] = 0
		pdmethSBoundary.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sBoundary','sEdge','id','int']
		pdmethSBoundary['Context'] = pdmethSBoundary['mstop'] + 1
		pdmethSBoundary['BackContext'] = pdmethSBoundary['mstart'] -1
		pdmethSBoundary['Nuc'] = simpleFasta(getFeatures(pdmethSBoundary[['mchr','BackContext','Context']].values.tolist()),faGenome)
		pdmethSBoundary['NucContext'] = simpleFasta(getFeatures(pdmethSBoundary[['mchr','mstop','Context']].values.tolist()),faGenome)
		pdmethSBoundary['NucCytosine'] = simpleFasta(getFeatures(pdmethSBoundary[['mchr','mstart','mstop']].values.tolist()),faGenome)
		pdmethSBoundary['NucBackContext'] = simpleFasta(getFeatures(pdmethSBoundary[['mchr','BackContext','mstart']].values.tolist()),faGenome)
		pdmethSBoundary['methLoc'] = pdmethSBoundary['int'].astype(int)+(pdmethSBoundary['mstart'].astype(int)-pdmethSBoundary['sBoundary'].astype(int))
		outpdmethSBoundary = pdmethSBoundary[['chr','mstart','mstop','sBoundary','sEdge','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outpdmethSBoundary.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
	else:
		outpdmethSBoundary = None

	methMiddle = methFeature.intersect(rangeFeatures[['chr','sCenter','eCenter','id']].values.tolist(),wb=True,wa=True)
	if len(methMiddle) != 0:
		pdmethFeature = bedtoolToPanda(methMiddle)
		pdmethFeature['int'] = (((num - uce)/2) + inuce)
		pdmethFeature.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sCenter','eCenter','id','int']
		pdmethFeature['Context'] = pdmethFeature['mstop'] + 1
		pdmethFeature['BackContext'] = pdmethFeature['mstart'] -1
		pdmethFeature['Nuc'] = simpleFasta(getFeatures(pdmethFeature[['mchr','BackContext','Context']].values.tolist()),faGenome)
		pdmethFeature['NucContext'] = simpleFasta(getFeatures(pdmethFeature[['mchr','mstop','Context']].values.tolist()),faGenome)
		pdmethFeature['NucCytosine'] = simpleFasta(getFeatures(pdmethFeature[['mchr','mstart','mstop']].values.tolist()),faGenome)
		pdmethFeature['NucBackContext'] = simpleFasta(getFeatures(pdmethFeature[['mchr','BackContext','mstart']].values.tolist()),faGenome)
		pdmethFeature['methLoc'] = pdmethFeature['int'].astype(int)+(pdmethFeature['mstart'].astype(int)-pdmethFeature['sCenter'].astype(int))
		outpdmethFeature = pdmethFeature[['chr','mstart','mstop','sCenter','eCenter','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outpdmethFeature.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
	else:
		outpdmethFeature = None

	methEBoundary = methFeature.intersect(rangeFeatures[['chr','eEdge','eBoundary','id']].values.tolist(),wb=True,wa=True)
	if len(methEBoundary) != 0:
		pdmethEBoundary = bedtoolToPanda(methEBoundary)
		pdmethEBoundary['int'] = num-1
		pdmethEBoundary.columns = ['mchr','mstart','mstop','methCov','methPer','chr','eEdge','eBoundary','id','int']
		pdmethEBoundary['Context'] = pdmethEBoundary['mstop'] + 1
		pdmethEBoundary['BackContext'] = pdmethEBoundary['mstart'] -1
		pdmethEBoundary['Nuc'] = simpleFasta(getFeatures(pdmethEBoundary[['mchr','BackContext','Context']].values.tolist()),faGenome)
		pdmethEBoundary['NucContext'] = simpleFasta(getFeatures(pdmethEBoundary[['mchr','mstop','Context']].values.tolist()),faGenome)
		pdmethEBoundary['NucCytosine'] = simpleFasta(getFeatures(pdmethEBoundary[['mchr','mstart','mstop']].values.tolist()),faGenome)
		pdmethEBoundary['NucBackContext'] = simpleFasta(getFeatures(pdmethEBoundary[['mchr','BackContext','mstart']].values.tolist()),faGenome)
		pdmethEBoundary['methLoc'] = pdmethEBoundary['int'].astype(int)-(pdmethEBoundary['eBoundary'].astype(int)-pdmethEBoundary['mstop'].astype(int))
		outpdmethEBoundary = pdmethEBoundary[['id','methPer','methLoc','methCov']]
		outpdmethEBoundary = pdmethEBoundary[['chr','mstart','mstop','eEdge','eBoundary','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outpdmethEBoundary.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
	else:
		outpdmethEBoundary = None

	methList = [outpdmethSBoundary,outpdmethFeature,outpdmethEBoundary]
	concatMeth = pd.concat(methList)
	sortMeth = concatMeth.sort_values(['methLoc'],ascending=True)
	return sortMeth

def compactMeth(mFiles,rangeFeatures,num,uce,inuce,methCovThresh,methPerThresh,faGenome):
	outMeth = []
	for methName in mFiles:
		methFeatures = eachFileProcess(methName)
		pdmethThresh = methThreshold(methFeatures,methCovThresh,methPerThresh)
		methPosition = methIntersect(rangeFeatures,pdmethThresh,num,uce,inuce,faGenome)
		methPosition['tissue'] = methName.replace('.bed','')
		stringDF = rangeFeatures[['id','combineString']]
		methMerge = pd.merge(methPosition,stringDF,how='left',on='id')
		methMerge['methLocBEnd'] = methMerge['methLoc'] - 1
		methMerge['methLocCEnd'] = methMerge['methLoc'] + 1
		methMerge['methLocEnd'] = methMerge['methLoc'] + 2
		methMerge['Cytosine'] = methMerge.apply(lambda row: row['combineString'][row['methLoc']:row['methLocCEnd']],axis=1)
		methMerge['Context'] = methMerge.apply(lambda row: row['combineString'][row['methLocCEnd']:row['methLocEnd']],axis=1)
		methMerge['BackContext'] = methMerge.apply(lambda row: row['combineString'][row['methLocBEnd']:row['methLoc']],axis=1)
		methMerge['ContextCheck'] = methMerge.apply(lambda row: row['combineString'][row['methLocBEnd']:row['methLocEnd']],axis=1)
		methMerge['methFreq'] = methMerge.groupby(['methLoc','Cytosine'])['methLoc'].transform('count')
		methMerge['Nuc'] = methMerge['Nuc'].str.upper()
		methMerge['sameSeq'] = methMerge['Nuc'] == methMerge['ContextCheck']
		
		# If the nucleotide in the cytosine column is 'G', make the context the other direction (reverse complement later, in graphing, in order to differentiate between strands)
		methMerge.loc[methMerge['Cytosine'] == 'G', 'Context'] = methMerge['BackContext']
		
		# sameSeq might be 'False' if 1) the c is at the end border for the downstream boundary, 2) the sequence bridges the sequence split for the upstream boundary
		falseMeth = (methMerge[methMerge['sameSeq'] == False])
		
		# Conditionally update contexts where False for matches between sequence and methylation nucleotides- c get context, and g gets backcontext
		methMerge.loc[methMerge['sameSeq'] == False,'Cytosine'] = methMerge['NucCytosine']
		methMerge.loc[(methMerge['sameSeq'] == False) & (methMerge['NucCytosine'] == 'C'),'Context'] = methMerge['NucContext']
		methMerge.loc[(methMerge['sameSeq'] == False) & (methMerge['NucCytosine'] == 'G'),'Context'] = methMerge['NucBackContext']
		
		print 'There are {0} instances at {1} where methylation context did not match between methylation bedfile and sequence in {2}'.format(len(falseMeth.index),falseMeth['methLoc'].tolist(),methName)
		subMeth = methMerge[['id','methLoc','methPer','methCov','methFreq','Cytosine','Context','tissue']]
		outMeth.append(subMeth)
	print 'Discrepancy of context are acceptable at the end of the sequence, or at the split between the middle and boundary, the context extracted from the original methylation file will be used instead'
	pdMeth = pd.concat(outMeth)
	return pdMeth

def main(mFiles,rangeFeatures,num,uce,inuce,methCovThresh,methPerThresh,faGenome):
	pdMeth = compactMeth(mFiles,rangeFeatures,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	return pdMeth

if __name__ == "__main__":
	main()