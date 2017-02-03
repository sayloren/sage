"""
    Script to create histograms for genomic elements based on AT content of bins
    
    Inputs:
    1 - txt file with the list of feature files to process, each bed file should have coordinates and type for fourth column
    2 - genome file with size of all chromosomes
    3 - fasta file for whole genome
    
    Outputs:
    1 - txt with the id of feature and the % AT content for each bin
    2 - histograms (both smooth, from Ruth's Plot_distances_IQR6 scirpt, and regular) for % AT content for each type of feature, and at each bin for each type of feature
    3 - stats for % AT content (modeled on Ruth's GC_summary_stats4.py)
    
    Wren Saylor
    January 13 2017
    
    Version 2
    Incoprorating Ruth's Plot_distance_IQR6.py, where bins are displayed continuously
    
    TODO:
    1 - line graphs total AT content at the edges of the feature ( think about the fang graphs )
    2 - try pybedtools alternative functions
    3 - make histograms for the sloped edges around the elements, slop
    4 - find way to keep id and coordinates
    5 - make 'all' print with the types on the same pdf for the boxplot
    
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

# 6 - stats by bin; use this file in Ruth's pipeline for smooth histograms
def smoothStats(pdFeatures, fileName):
    noIDFeatures = pdFeatures.drop(pdFeatures.columns[[0]],axis=1)
        transposeFeatures = noIDFeatures.transpose()
        transFeatures = transposeFeatures.astype(int)
        transFeatures['0'] = transFeatures.mean(axis=1)
        transFeatures['lower_quartile'] = transFeatures.quantile(axis=1,q=0.25)
        transFeatures['median'] = transFeatures.median(axis=1)
        transFeatures['upper_quartile'] = transFeatures.quantile(axis=1,q=0.75)
        outTrans = transFeatures[['0','upper_quartile','median','lower_quartile']]
        outTrans.to_csv('SmoothStatsfor_all_{0}.txt'.format(fileName),index=True,sep="\t")
        smoothHistplot(outTrans, 'Distances_with_randomIQR_all_{0}'.format(fileName))
        
        # By type
        typeList = ['exonic','intronic','intergenic']
        for element in typeList:
            boolType = pdFeatures[pdFeatures[0] == element]
                transType = boolType.loc[:,1:].transpose()
                intType = transType.astype(int)
                intType['0'] = intType.mean(axis=1)
                intType['lower_quantile'] = intType.quantile(axis=1,q=0.25)
                intType['upper_quantile'] = intType.quantile(axis=1,q=0.75)
                intType['median'] = intType.median(axis=1)
                outType = intType[['0','lower_quantile','median','upper_quantile']]
                outType.to_csv('SmoothStatsfor_{0}_{1}.txt'.format(element,fileName),index=True,sep="\t")
                smoothHistplot(outType, 'Distances_with_randomIQR_{0}_{1}'.format(element,fileName))

# 6 - make regular histogram
def makeHistograms(name, fileName, group):
    strFilename = 'Hist_full_result_for_{0}_{1}.pdf'.format(name,fileName)
        plt.figure()
        g = group.hist() # by = type
        plt.savefig(strFilename, format='pdf', bbox_inches='tight')
        plt.clf()

# 7 - make boxplots
def makeBoxplots(fileName, group):
    strFilename = 'Boxplot_full_result_for_{0}.pdf'.format(fileName)
        plt.figure()
        g = group.groupby([0]).boxplot(subplots=True,return_type='axes')
        plt.savefig(strFilename,format='pdf', bbox_inches='tight')
        plt.clf()

# # 8 - Ruth's smooth histogram formating for labels and label locations
# def getMiddleBinsGetLabels(pdValues):
# 	seriesXValues = pdValues.ix[:, 0]
# 	Get distance between bin
# 	intFirstBinStart = seriesXValues.iloc[0]
# 	intSecondBinStart = seriesXValues.iloc[1]
# 	intInterBinDistance = intSecondBinStart - intFirstBinStart
# 	intHalfBinDistance = int(intInterBinDistance/2)
# 	Add intHalfBinDistance to each x value
# 	pdValues.columns = ['bin_start','lower_quartile','median','upper_quartile']
# 	pdValues['mid_bin'] = pdValues.apply(lambda row: row['bin_start'] + intHalfBinDistance, axis=1)
# 	pdWithMidBins = pdValues[['mid_bin','lower_quartile','median','upper_quartile']]
# 	print pdWithMidBins
# 	Create list of axis labels
# 	arLabels = []
# 	arMidBins = pdWithMidBins['mid_bin'].values.tolist()
# 	for floatMidBin in arMidBins:
# 		intLower = (floatMidBin-intHalfBinDistance)/1000
# 		intUpper = (floatMidBin + intHalfBinDistance)/1000
# 		strLabel = '{0} to {1}'.format(intLower, intUpper)
# 		arLabels.append(strLabel)
# 	Thin out the labels if necessary
# 	intNumberLabels = len(arLabels)
# 	arThinnedLabels = []
# 	arLabelPositions = []
# 	intThin = 1
# 	for n in range(0, intNumberLabels, intThin):
# 		thinLabel = arLabels[n]
# 		arThinnedLabels.append(thinLabel)
# 		floatLabelPosition = arMidBins[n]
# 		arLabelPositions.append(floatLabelPosition)
# 	print arLabelPositions
# 	return pdWithMidBins, arThinnedLabels, arLabelPositions

# 9 - Ruth's smooth histogram plotting
def smoothHistplot(pdValues, strName):
    npX = range(0,len(pdValues))
        npQuery = pdValues.ix[:, 0].values
        npUpperQ = pdValues.ix[:, 1].values
        npMedian = pdValues.ix[:, 2].values
        npLowerQ = pdValues.ix[:, 3].values
        tupXlim = (0,len(pdValues))
        intXLower = tupXlim[0]
        intXUpper = tupXlim[1]
        plt.figure(figsize=(6, 4))
        sns.set_style('ticks')
        ax = plt.plot(npX, npQuery, linewidth=1, color='orange', label='')
        plt.plot(npX, npQuery, marker='o', markersize=7, color='orange', clip_on=False)
        plt.fill_between(npX, npUpperQ, npLowerQ, facecolor='#bdbdbd', edgecolor='#bdbdbd', label = '')
        sns.despine(offset=5)
        #sns.axlabel('Distance to nearest breakpoint (Mb)', 'Frequency', fontsize=10)
        plt.xticks(rotation=45, ha='right')
        plt.xlim(intXLower, intXUpper)
        plt.savefig('{0}.pdf'.format(strName), format='pdf', bbox_inches='tight')
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
                
                # Make smooth histogram
                smoothStats(pdFeatures, fileName)
                
                # Make files for AT stats
                arStatColumns = ['NumObs', 'floatMean', 'float25thPercentile', 'float75thPercentile', 'floatMin', 'floatMax']
                pdAllStats = pd.DataFrame(data=arArStats, index=indexStats, columns=arStatColumns)
                pdAllStats.to_csv('ATStatsfor_{0}.txt'.format(fileName), sep='\t')
                
                            makeBoxplots(fileName, pdFeatures)

if __name__ == "__main__":
    main()
