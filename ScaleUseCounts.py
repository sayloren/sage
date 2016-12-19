"""

Script for taking in a bed file with the number of overlaps of regions of interest as 
columns, counting the number of overlaps for each column in regions of interest, 
then scaling that range for the useScore option in the UCSC genome browser.

Input: A bedfile, with descriptive headers of form "Chromosome \t Start \t Stop \t ID, and other columns

Output:

To Do:
-Allow so that desired columns to sum over are named in the argparse input
-Also, header is printed, to view with UCSC, will need to have the useScale = 1 option,
and whatever other descriptions user wants added
-Possibly use pybedtools.featurefuncs.bedgraph_scale

Wren Saylor
December 2 2016

"""

import argparse
import pybedtools as pbt
import pandas as pd


def get_args():
	parser = argparse.ArgumentParser(description='Script for taking in a bed file with the '
	'number of overlaps of regions of interest as columns, counting the number of overlaps '
	'for each column in regions of interest, then scaling that range for the useScore option '
	'in the UCSC genome browser.')
	parser.add_argument("file", type=argparse.FileType('rU'))
	return parser.parse_args()

# 1 - read in the file
def getFeatures(fileName):
	btFeatures = pbt.BedTool(fileName)
	return btFeatures
	
# 2 - turn into panda
def bedtoolToPanda(btFeatures):
	pdObject = pd.read_table(btFeatures.fn, header=0)
	return pdObject

# 3 - count the overlaps
def pandaColSums(pdObject):
	pdObject['TotalCount'] = pdObject.iloc[:,4::].sum(axis=1) # if matrix is already summed, comment this line out
	pdObjectCounts = pdObject[['Chromosome','Start','Stop','ID','TotalCount']]
	return pdObjectCounts

# 4 - scale range 0-945, from UCSC's useScale=1 bed file option
# Formula for scaling: NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
def scaleTotalCount(pdObjectCounts):
	pdOMin = pdObjectCounts['TotalCount'].min()
	pdOMax = pdObjectCounts['TotalCount'].max()
	pdObjectCounts['ScaledCount'] = (((pdObjectCounts['TotalCount'] - pdOMin) * (495 - 0)) / (pdOMax - pdOMin) + 0).astype(int)
	pdObjScaled = pdObjectCounts[['Chromosome','Start','Stop','ID','ScaledCount']]
	print pdObjScaled
	return pdObjScaled

# save panda file
def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', index=False)

# turn into bedtool (didn't use)
def pandaToBedtool(panda): # columns 5+
	arArFeatures = panda.values.tolist()
	btoutFeatures = getFeatures(arArFeatures)
	return btoutFeatures

# save bedtool file (didn't use)
def saveBedTool(btoutFeatures, strFilename):
	btoutFeatures.saveas(strFilename)


def main(args):

	aFiles = [line.strip() for line in args.file]
	for fileName in aFiles:
		btFeatures = getFeatures(fileName)
		pdObject = bedtoolToPanda(btFeatures)
		pdObjectCounts = pandaColSums(pdObject)
		pdObjScaled = scaleTotalCount(pdObjectCounts)
		savePanda(pdObjScaled,"scaledTotaled_" + str(fileName))
		#btoutFeatures = pandaToBedtool(pdObjScaled) # (didn't use)
		#saveBedTool(btoutFeatures,"scaledTotaled_" + str(fileName)) # (didn't use)
		

if __name__ == "__main__":
	args = get_args()
	main(args)
