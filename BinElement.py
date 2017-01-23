"""
Script to create histograms for genomic elements based on AT content of bins

Inputs:
1 - txt file with the list of feature files to process, each bed file should have coordinates and type for fourth column
2 - genome file with size of all chromosomes
3 - fasta file for whole genome

Outputs: 
1 - txt with the id of feature and the % AT content for each bin
2 - histogram for % AT content for each type of feature
3 - histogram for % AT content at each bin for each type of feature
4 - stats for % AT content ( modeled on Ruth's GC_summary_stats4.py )

Wren Saylor
January 13 2017

TODO:
1 - line graphs total AT content at the edges of the feature ( think about the fang graphs )
2 - try pybedtools alternative functions
3 - make histograms for the sloped edges around the elements, slop
4 - find way to keep id and coordinates
5 - make 'all' print with the types on the same pdf
"""

import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pybedtools as pbt
import Bio
from Bio import SeqIO
from Bio import Seq

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a list of paths to the files separated by newlines')
	parser.add_argument('-c', '--chunks', metavar = 'chunks', required=True, type = int, help = "Number of chunks to break up the UCE into")
# GENOME FILES
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
	parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")
	return parser.parse_args()

# 1 - read in files
def eachFileProcess(fileName):
	btFeatures = pbt.BedTool(fileName)
	return btFeatures
	
# 2 - get sequence
def getFasta(btFeatures,faGenome,fileName):
	saveFeatures = 'Seq_type_result_for_{0}.txt'.format(fileName)
	seqFeatures = btFeatures.sequence(fi=faGenome,fo=saveFeatures,name=True) # to get the coordinates, remove name=True, wont pass through script, though
	return saveFeatures

# 3 - chunk into bins
def getChunks(seq_record, sizeChunk):
	currentID, currentFeature = seq_record.id[:], seq_record.seq[:]
	chunkFeatures = []
	chunkFeatures.append(currentID)
	n = len(currentFeature)
	avg = n / sizeChunk
	remainders = n % sizeChunk
	start, end = 0, avg
	while start < n:
		if remainders > 0:
			end = end + 1
			remainders = remainders -1
		currentChunk = currentFeature[start:end]
		percentage = eval('100 * float(currentChunk.count("A") + currentChunk.count("T") + currentChunk.count("a") + currentChunk.count("t")) / len(currentChunk)')
		chunkFeatures.append(percentage)
		start, end = end, end + avg
	return chunkFeatures # avg

# 4a - save file from panda
def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', header=False, index=False)

# 4b - save file from bedtool
def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

# TO DO (alt)- chunk into bins with pbt
# def pyGetChunks(btFeatures, sizeChunk, faGenome):
# 	winFeatures = pbt.BedTool().window_maker(b=str(btFeatures), w = sizeChunk)
# 	faFeatures = getFasta(winFeatures,faGenome)
# 	btNucResult = faFeatures.nucleotide_content(fi=faGenome)
# 	return btNucResult

# TO DO (alt) - get sequences with pbt
# def getpdGC(btNucResult):
# 	#Make bedtool into pandas object
# 	pandasNuc = pd.read_table(btNucResult.fn)
# 	#Get the name of the column that holds the information for GC fraction
# 	strGCColumn = getGCColumnName(pandasNuc)
# 	#Get a pandas dataframe containing just the GC fraction
# 	pdGC = pandasNuc[strGCColumn]
# 	return pdGC

# 5 - stats file; use this stats data to create a position graph for the UCEs
def getStats(collectFeatures):
	for column in collectFeatures:
		intN = collectFeatures[column].count()
		floatMean = collectFeatures[column].mean()
		float25thPercentile = collectFeatures[column].quantile(q=0.25)
		float75thPercentile = collectFeatures[column].quantile(q=0.75)
		floatMin = collectFeatures[column].min()
		floatMax = collectFeatures[column].max()
		arStats = [intN, floatMean, float25thPercentile, float75thPercentile, floatMin, floatMax]
		return arStats

# 6 - make histogram
def makeHistograms(name, fileName, group):
	strFilename = 'Hist_full_result_for_{0}_{1}.pdf'.format(name,fileName)
	plt.figure()
	g = group.hist() # by = type
	plt.savefig(strFilename, format='pdf', bbox_inches='tight')
	plt.clf()
	
def makeBoxplots(fileName, group):
	strFilename = 'Boxplot_full_result_for_{0}.pdf'.format(fileName)
	plt.figure()
	g = group.groupby([0]).boxplot(subplots=True,return_type='axes')
	plt.savefig(strFilename,format='pdf', bbox_inches='tight')
	plt.clf()

def main():
	args = get_args()
	aFiles = [line.strip() for line in args.file]
	sizeChunk = args.chunks
	sizeGenome = args.genome
	faGenome = args.fasta
	# for file in master txt file
	for fileName in aFiles:
		btFeatures = eachFileProcess(fileName)
		saveFeatures = getFasta(btFeatures,faGenome,fileName)
		allFeatures = []
		for seq_record in SeqIO.parse(saveFeatures, "fasta"):
			chunkFeatures = getChunks(seq_record, sizeChunk)
			allFeatures.extend([chunkFeatures])
		pdFeatures = pd.DataFrame(allFeatures)
		savePanda(pdFeatures,"Bin_type_result_for_{0}.txt".format(fileName))
	
		# All UCEs, by bin
		collectFeatures = pdFeatures.drop(pdFeatures.columns[0],axis=1)
		stackFeatures = collectFeatures.stack()
		arArStats = []
		indexStats = ['all']
		allStats = getStats(collectFeatures)
		arArStats.append(allStats)
		
		# Each type of UCE, by bin
		valFeatures = pdFeatures.groupby([0])
		for name, group in valFeatures:
			group = group.drop(group.columns[0],axis=1)
			stack = group.stack()
			makeHistograms(name, fileName, stack)
			makeHistograms(name+"_byBin", fileName, group)
			typeStats = getStats(group)
			arArStats.append(typeStats)
			indexStats.append(name)
			
		arStatColumns = ['NumObs', 'floatMean', 'float25thPercentile', 'float75thPercentile', 'floatMin', 'floatMax']
		pdAllStats = pd.DataFrame(data=arArStats, index=indexStats, columns=arStatColumns)
		pdAllStats.to_csv('ATStatsfor_{0}.txt'.format(fileName), sep='\t')
		makeBoxplots(fileName, pdFeatures)

if __name__ == "__main__":
	main()
