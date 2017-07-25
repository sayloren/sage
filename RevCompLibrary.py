"""
Script to perform RC sorting

Wren Saylor
July 5 2017
"""
import argparse
import pandas as pd
from FangsLibrary import compactWindow
from ElementLibrary import eachFileProcess
from MethylationLibrary import compactMeth

# Methylation RCsorting
def methDirection(negStr,posStr,mFiles,num,uce,inuce,methCovThresh,methPerThresh,faGenome):

	posMeth = compactMeth(mFiles,posStr,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	negMeth = compactMeth(mFiles,negStr,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	
	# Zip reversed range to make a dictionary for replacing the location of the neg methylation
	originalRange = range(0,num)
	reverseRange = originalRange[::-1]
	rangeDict = dict(zip(originalRange,reverseRange))
	
	# Zip reverse complement sequence for replacing the nucleotides for neg methylation
	seqDict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

	# Convert neg Meth df
	negMeth['methLocNew'] = negMeth.methLoc.map(rangeDict)
	negMeth['CytosineNew'] = negMeth.Cytosine.map(seqDict)
	negMeth['ContextNew'] = negMeth.Context.map(seqDict)
	negMethNew = negMeth[['id','methLocNew','methPer','methCov','methFreq','CytosineNew','ContextNew','tissue']]
	negMethNew.columns = ['id','methLoc','methPer','methCov','methFreq','Cytosine','Context','tissue']
	
	# Concat pos and revised neg meth dfs
	frames = [posMeth,negMethNew]
	catMerge = pd.concat(frames)
	
	# Update Frequencey count column
	catMerge['methFreqNew'] = catMerge.groupby(['methLoc','tissue','Cytosine'])['methLoc'].transform('count')
	outMerge = catMerge[['id','methLoc','methPer','methCov','methFreqNew','Cytosine','Context','tissue']]
	outMerge.columns = ['id','methLoc','methPer','methCov','methFreq','Cytosine','Context','tissue']

	return outMerge

# Sliding window RCsorting
def slideDirection(negStr,posStr,num,uce,inuce,window,nucLine):
	negDF, negNames = compactWindow(negStr['reverseComplement'],negStr['id'],num,uce,inuce,window,nucLine)
	posDF, posNames = compactWindow(posStr['combineString'],posStr['id'],num,uce,inuce,window,nucLine)
	compWindow = []
	for x, y in zip(negDF, posDF):
		tempCat = pd.concat([x,y],axis=1)
		tempGroup = tempCat.groupby(tempCat.columns,axis=1).sum()
		compWindow.append(tempGroup)
	return compWindow, negNames # or could be posNames, they will be the same

# Separate on plus and minus orientation, RCsort and return methylation and sliding window computations
def dirLine(directionFeatures,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):
	negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
	posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
	compWindow, compNames = slideDirection(negStr,posStr,num,uce,inuce,window,nucLine)
	if any(x in graphs for x in ['methylation','cluster']):
		groupMeth = methDirection(negStr,posStr,mFiles,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(directionFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):
	groupMeth,compWindow,compNames = dirLine(directionFeatures,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
	print 'Completed reverse complement sorting for {0} items, with {1} bin sorting'.format(len(directionFeatures.index),binDir)
	return groupMeth,compWindow,compNames

if __name__ == "__main__":
	main()