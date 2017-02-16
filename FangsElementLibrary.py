"""
    Script to make fangs and do some further proccessing
    
    Inputs:
    1 - txt file with the list of feature files to process, each bed file should have coordinates and type for fourth column
    2 - genome file with size of all chromosomes
    3 - fasta file for whole genome
    
    Outputs:
    Fasta for features
    Graphs for Fangs
    File for directionality
    
    Wren Saylor
    January 13 2017
    
    To Do:
    Make the passing to the graphs more neat
    Print number of - + = on plot
    Separate A T
    Stats for significance
    Global organization
    Middle data
    Diff each UCE, graph, then graph stat for diffs of each type
    Use rate of change to get fang threshold for each UCE, return number of UCEs with 1/2/0 fangs
    label graphs (fang location, significance)
    try chaning sliding window size
    
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
import re
import numdifftools as nd
import matplotlib.gridspec as gridspec
# from matplotlib_scalebar.scalebar import ScaleBar

def get_args():
    parser = argparse.ArgumentParser(description="Description")
        parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a list of paths to the files separated by newlines')
        # GENOME FILES
        parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
        parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")
        return parser.parse_args()

# 1 - read in files
def eachFileProcess(fileName):
    btFeatures = pbt.BedTool(fileName)
        return btFeatures

# 2 - get bt features
def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
        return btFeatures

# 3 - get the correct range for fang evaluation
def getRange(btFeatures,fileName):
    # for each type in list of element type
    midFeatures = pd.read_table(btFeatures.fn, header=None)
        midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
        midFeatures['sCenter'] = midFeatures['middle'] - 50
        midFeatures['eCenter'] = midFeatures['middle'] + 50
        midFeatures['sEdge'] = midFeatures.loc[:,1] + 50
        midFeatures['eEdge'] = midFeatures.loc[:,2] - 50
        midFeatures['sBoundary'] = midFeatures.loc[:,1] - 200
        midFeatures['eBoundary'] = midFeatures.loc[:,2] + 200
        midFeatures['start'] = midFeatures.loc[:,1]
        midFeatures['end'] = midFeatures.loc[:,2]
        midFeatures['type'] = midFeatures.loc[:,3]
        midFeatures['chr'] = midFeatures.loc[:,0]
        rangeFeatures = midFeatures[['type','chr','sBoundary', 'start', 'sEdge', 'sCenter','eCenter','eEdge','end','eBoundary']]
        return rangeFeatures

# 4 - get fasta strings for each desired region
def getFasta(btFeatures,faGenome,fileName):
    # Just the fasta files for the UCEs
    saveAll = 'Seq_result_for_all_{0}.txt'.format(fileName)
        seqAll = btFeatures.sequence(fi=faGenome,fo=saveAll)
        saveExonic = 'Seq_results_for_exonic_{0}.txt'.format(fileName)
        seqExonic = btFeatures.sequence(fi=faGenome,fo=saveExonic).filter(lambda x: x[name] == 'exonic')
        saveIntronic = 'Seq_results_for_intronic_{0}.txt'.format(fileName)
        seqIntronic = btFeatures.sequence(fi=faGenome,fo=saveIntronic).filter(lambda x: x[name] == 'intronic')
        saveIntergenic = 'Seq_results_for_intergenic_{0}.txt'.format(fileName)
        seqIntergenic = btFeatures.sequence(fi=faGenome,fo=saveIntergenic).filter(lambda x: x[name] == 'intergenic')
        return saveAll, saveExonic, saveIntronic, saveIntergenic

# 5 - used in btRange to extract just the fasta strings
def simpleFasta(inFeature,faGenome):
    seqFeature = inFeature.sequence(fi=faGenome)
        outFeature = pd.read_table(seqFeature.seqfn)
        outSequence = outFeature[::2]
        outSequence = outSequence.reset_index(drop=True)
        return outSequence

# 6 - get the strings for sliding window regions
def btRange(rangeFeatures,faGenome):
    # Now the fasta files with the ranged regions
    rangeFeatures['sBoundarySeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','sBoundary','start']].values.tolist()),faGenome)
        rangeFeatures['sEdgeSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','start','sEdge']].values.tolist()),faGenome)
        rangeFeatures['MiddleSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','sCenter','eCenter']].values.tolist()),faGenome)
        rangeFeatures['eEdgeSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','eEdge','end',]].values.tolist()),faGenome)
        rangeFeatures['eBoundarySeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','end','eBoundary']].values.tolist()),faGenome)
        rangeFeatures['feature'] = simpleFasta(getFeatures(rangeFeatures[['chr','start','end']].values.tolist()),faGenome)
        #rangeFeatures['combineString'] = rangeFeatures[['sBoundarySeq','sEdgeSeq','MiddleSeq','eEdgeSeq','sBoundarySeq']].apply(lambda x: ''.join(x),axis=1)
        rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
        rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
        return rangeFeatures

# 7 - sliding window
def slidingWindow(element):
    winFeatures = []
        winA = []
        winT = []
        winG = []
        winC = []
        n = 600 # take out hard coding
        start, end = 0, 11 # possibly make larger sliding window?
        while end < n:
            current = element[start:end]
                percentageAT = eval('100 * float(current.count("A") + current.count("T"))/ len(current)')
                perA = eval('100 * float(current.count("A"))/ len(current)')
                perT = eval('100 * float(current.count("T"))/ len(current)')
                perG = eval('100 * float(current.count("G"))/ len(current)')
                perC = eval('100 * float(current.count("C"))/ len(current)')
                winFeatures.append(percentageAT)
                winA.append(perA)
                winT.append(perT)
                winG.append(perG)
                winC.append(perC)
                start, end = start + 1, end + 1
        # print start, end  # a check, to make sure that the string is the correct size
    return winFeatures, winA, winT, winG, winC

# 8 - make some graphs!
def endLinegraphs(pdWindow,pdA,pdT,pdG,pdC,fileName):
    fillX = range(0,589) # not 600 because of the sliding window
        sns.set_style('ticks')
        fig = plt.figsize=(14,7)
        gs = gridspec.GridSpec(5,1,height_ratios=[2,3,3,1,1])
        gs.update(hspace=.5)
        info = str(fileName) + ', '+ str(len(pdWindow)) + ' - ' "UCES"
        plt.suptitle(info,fontsize=14)
        
        # Plot the mean AT content with a std of 1
        ax2 = plt.subplot(gs[1])
        #ax2.update(hspace=1)
        ax2.plot(fillX,pdWindow.mean(axis=0),linewidth=1, color='#081d58')
        ax2.fill_between(fillX,pdWindow.mean(axis=0)+pdWindow.std(axis=0,ddof=1), pdWindow.mean(axis=0)-pdWindow.std(axis=0,ddof=1),facecolor = '#081d58',label='',alpha=0.4)
        ax2.axvline(x=195,linewidth=5,color='#ae017e',alpha=0.3)
        ax2.axvline(x=395,linewidth=5,color='#ae017e',alpha=0.3)
        ax2.axvline(x=195,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax2.axvline(x=395,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax2.set_xticks([])
        ax2.set_yticks([30,60,90])
        ax2.set_title('Mean AT Content',size=10)
        # 	scalebar = ScaleBar(0.01) # 1 pixel = 0.2 meter
        # 	plt.gca().add_artist(scalebar)
        
        
        # Plot the std = 1
        ax1 = plt.subplot(gs[0])
        #ax1.update(hspace=1)
        ax1.plot(fillX,pdWindow.std(axis=0,ddof=1),linewidth=1, color='#081d58')
        ax1.axvline(x=195,linewidth=5,color='#ae017e',alpha=0.3)
        ax1.axvline(x=395,linewidth=5,color='#ae017e',alpha=0.3) # not sure how to set the size to that of the sliding window
        ax1.axvline(x=195,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax1.axvline(x=395,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax1.set_xticks([])
        ax1.set_yticks([14,17,20])
        ax1.set_title('Standard Deviation, DOF=1',size=10)
        
        # Nucleotides
        ax3 = plt.subplot(gs[2])
        ax3.plot(fillX,pdA.mean(axis=0),linewidth=1, color='#810f7c',label='A') # purple
        ax3.plot(fillX,pdT.mean(axis=0),linewidth=1, color='#006d2c',label='T') # green
        ax3.plot(fillX,pdG.mean(axis=0),linewidth=1, color='#0868ac',label='G') # blue
        ax3.plot(fillX,pdC.mean(axis=0),linewidth=1, color='#b30000',label='C') # red
        ax3.axvline(x=195,linewidth=5,color='#ae017e',alpha=0.3)
        ax3.axvline(x=395,linewidth=5,color='#ae017e',alpha=0.3) # not sure how to set the size to that of the sliding window
        ax3.axvline(x=195,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax3.axvline(x=395,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax3.set_xticks([])
        ax3.set_yticks([10,25,40])
        ax3.set_title('Mean Nucleotide Content Separately',size=10)
        #ax3.legend(loc=0,fontsize=5)
        
        # Plot the median and quantiles
        # 	ax3.plot(fillX,pdWindow.median(axis=0),linewidth=1, color='#081d58')
        # 	ax3.fill_between(fillX,pdWindow.min(axis=0),pdWindow.max(axis=0),facecolor='#edf8b1', label='')
        # 	ax3.fill_between(fillX,pdWindow.quantile(axis=0,q=.05),pdWindow.quantile(axis=0,q=0.95),facecolor = '#c7e9b4',label='')
        # 	ax3.fill_between(fillX,pdWindow.quantile(axis=0,q=0.25),pdWindow.quantile(axis=0,q=0.75),facecolor = '#7fcdbb',label='')
        # 	ax3.axvline(x=195,linewidth=5,color='#ae017e',alpha=0.5)
        # 	ax3.axvline(x=395,linewidth=5,color='#ae017e',alpha=0.5)
        
        # Plot the derivative
        ax4 = plt.subplot(gs[3])
        #ax3.updatae(hspace=1)
        ax4.plot(fillX,pdWindow.mean(axis=0).diff(periods=11),linewidth=1, color='#081d58')
        ax4.axvline(x=195,linewidth=5,color='#ae017e',alpha=0.3)
        ax4.axvline(x=395,linewidth=5,color='#ae017e',alpha=0.3)
        ax4.axvline(x=195,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax4.axvline(x=395,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax4.set_xticks([])
        ax4.set_yticks([-10,0,10])
        ax4.set_title('AT Content Rate of Change',size=10)
        
        # Plot the second derivative
        # Convolution kernel equal to second derivative
        #http://stackoverflow.com/questions/13691775/python-pinpointing-the-linear-part-of-a-slope
        #http://stackoverflow.com/questions/9300430/algorithm-is-there-a-linear-trend-in-data?rq=1
        smooth = 50
        x1 = np.linspace(-3,3,smooth)
        y1 = (4*x1**2-2) * np.exp(-x1**2) / smooth * 8
        mean = pdWindow.mean(axis=0)
        y_conv = np.convolve(mean,y1,mode="same")
        ax5 = plt.subplot(gs[4])
        #ax4.update(hspace=1)
        ax5.plot(fillX,y_conv,linewidth=1, color='#081d58')
        ax5.set_ylim(-11,11)
        ax5.axvline(x=195,linewidth=5,color='#ae017e',alpha=0.3)
        ax5.axvline(x=395,linewidth=5,color='#ae017e',alpha=0.3)
        ax5.axvline(x=195,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax5.axvline(x=395,linewidth=.5,linestyle='dashed',color='#ae017e')
        ax5.set_yticks([-14,0,14])
        ax5.set_title('AT Content Second Derivative',size=10)
        ax5.set_xlabel('Base Pairs',size=14)
        # this is off for some reason????
        
        sns.despine(offset=3)
        plt.xlim(0,589)
        plt.savefig('Fangs_{0}.pdf'.format(fileName), format='pdf', bbox_inches='tight')
        plt.clf()

def dirLine(rangeFeatures,fileName):
    #posStr = (rangeFeatures[(rangeFeatures['compareFive'] == "+") & (rangeFeatures['compareTen'] == "+") & (rangeFeatures['compareTwenty'] == "+") & (rangeFeatures['compareThirty'] == "+") & (rangeFeatures['compareFourty'] == "+")])
    posStr = (rangeFeatures[(rangeFeatures['compareThirty'] == "+")])
        # 	posStr = (rangeFeatures[(rangeFeatures['compareTen'] == "+") &
        # 	(rangeFeatures['compareTwenty'] == "+") &
        # 	(rangeFeatures['compareThirty'] == "+")])
        outPos = []
        for element in posStr['combineString']:
            posFeature = slidingWindow(element)
                outPos.append(posFeature)
        posWindow = pd.DataFrame(outPos)
        endLinegraphs(posWindow,'Pos_30_{0}'.format(fileName))
        
        #negStr = (rangeFeatures[(rangeFeatures['compareFive'] == "-") & (rangeFeatures['compareTen'] == "-") & (rangeFeatures['compareTwenty'] == "-") & (rangeFeatures['compareThirty'] == "-") & (rangeFeatures['compareFourty'] == "-")])
        negStr = (rangeFeatures[(rangeFeatures['compareThirty'] == "-")])
        # 	negStr = (rangeFeatures[(rangeFeatures['compareTen'] == "-") &
        # 	(rangeFeatures['compareTwenty'] == "-") &
        # 	(rangeFeatures['compareThirty'] == "-")])
        outNeg = []
        #negStr = (rangeFeatures[rangeFeatures['compareFive'] == '-'])
        for element in negStr['combineString']:
            negFeature = slidingWindow(element)
                outNeg.append(negFeature)
        negWindow = pd.DataFrame(outNeg)
        endLinegraphs(negWindow,'Neg_30_{0}'.format(fileName))

#equStr = (rangeFeatures[(rangeFeatures['compareFive'] == "+") & (rangeFeatures['compareTen'] == "+") & (rangeFeatures['compareTwenty'] == "+") & (rangeFeatures['compareThirty'] == "+") & (rangeFeatures['compareFourty'] == "+") | (rangeFeatures['compareFive'] == "-") & (rangeFeatures['compareTen'] == "-") & (rangeFeatures['compareTwenty'] == "-") & (rangeFeatures['compareThirty'] == "-") & (rangeFeatures['compareFourty'] == "-")])
# 	dupStr = (rangeFeatures[((rangeFeatures['compareTen'] == "+") &
# 	(rangeFeatures['compareTwenty'] == "+") &
# 	(rangeFeatures['compareThirty'] == "+")) |
# 	((rangeFeatures['compareTen'] == "-") &
# 	(rangeFeatures['compareTwenty'] == "-") &
# 	(rangeFeatures['compareThirty'] == "-"))])
dupStr = (rangeFeatures[(rangeFeatures['compareThirty'] == "+") | (rangeFeatures['compareThirty'] == "-")])
    outDup = []
        #equStr = (rangeFeatures[rangeFeatures['compareFive'] == '='])
        for element in dupStr['combineString']:
            dupFeature = slidingWindow(element)
                outDup.append(dupFeature)
        dupWindow = pd.DataFrame(outDup)
        endLinegraphs(dupWindow,'reComb_30_{0}'.format(fileName))
        
        equStr = (rangeFeatures[(rangeFeatures['compareThirty'] == "=")])
        # 	(rangeFeatures['compareTwenty'] == "=") &
        # 	(rangeFeatures['compareThirty'] == "=")])
        outEqu = []
        for element in equStr['combineString']:
            equFeature = slidingWindow(element)
                outEqu.append(equFeature)
        equWindow = pd.DataFrame(outEqu)
                    endLinegraphs(equWindow,'Equ_30_{0}'.format(fileName)) # ,outA,outT

# 9 - save file from panda
def savePanda(pdData, strFilename):
    pdData.to_csv(strFilename, sep='\t', header=False, index=False)

# 9 - out put directionality, as inferred by comparing first and last n base pairs
def compareN(element,size):
    start = element[:size]
        end = element[-size:]
        perSize = []
        perSize.append(eval('100*int(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
        perSize.append(eval('100*int(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
        #outList.append(size) # uncomment if want the size of the end in the dataframe
        #outList.append(perSize) # uncomment if want the values compared in the dataframe
        if perSize[0] > perSize[1]: outList = '+'
        if perSize[1] > perSize[0]: outList = '-'
        if perSize[1] == perSize[0]: outList = '='
        return outList

def main():
    
    # Collect arguments
    args = get_args()
        aFiles = [line.strip() for line in args.file]
        sizeGenome = args.genome
        faGenome = args.fasta
        
        # Collect and processes files
        for fileName in aFiles:
            btFeatures = eachFileProcess(fileName)
                rangeFeatures = getRange(btFeatures, fileName)
                saveAll, saveExonic, saveIntronic, saveIntergenic = getFasta(btFeatures,faGenome,fileName)
                rangeFeatures = btRange(rangeFeatures,faGenome)
                outWindow = []
                outA = []
                outT = []
                outG = []
                outC = []
                for element in rangeFeatures['combineString']:
                    outFeature, winA, winT, winG, winC = slidingWindow(element)
                        outWindow.append(outFeature)
                        outA.append(winA)
                        outT.append(winT)
                        outG.append(winG)
                        outC.append(winC)
                pdWindow = pd.DataFrame(outWindow)
                pdA = pd.DataFrame(outA)
                pdT = pd.DataFrame(outT)
                pdG = pd.DataFrame(outG)
                pdC = pd.DataFrame(outC)
                endLinegraphs(pdWindow,pdA,pdT,pdG,pdC,fileName)
                
                rangeFeatures['compareFive'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],5)),axis=1)
                rangeFeatures['compareFive'] = rangeFeatures['compareFive'].astype(str)
                rangeFeatures['compareTen'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],10)),axis=1)
                rangeFeatures['compareTwenty'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],20)),axis=1)
                rangeFeatures['compareThirty'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],30)),axis=1)
                rangeFeatures['compareFourty'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],40)),axis=1)
                compareEnds = pd.DataFrame(rangeFeatures[['chr','start','end','compareFive','compareTen','compareTwenty','compareThirty','compareFourty']])
            savePanda(compareEnds,'Directionality_results_for_all_{0}.txt'.format(fileName))

#dirLine(rangeFeatures,fileName)

# 		typeList = ['exonic','intronic','intergenic']
# 		for type in typeList:
# 			boolType = (rangeFeatures[rangeFeatures['type'] == type])
# 			outType = []
# 			for element in boolType['combineString']:
# 				outFeature = slidingWindow(element)
# 				outType.append(outFeature)
# 			pdWindow = pd.DataFrame(outType)
# 			endLinegraphs(pdWindow,pdA,pdT,'{0}_{1}'.format(type,fileName))
#dirLine(rangeFeatures,'{0}_{1}'.format(type,fileName))


if __name__ == "__main__":
    main()
