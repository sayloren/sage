"""
    Script to make fangs and do some further proccessing
    
    Inputs:
    1 - txt file with the list of feature files to process, each bed file should have coordinates, id for fourth and type for fifth column
    2 - genome file with size of all chromosomes
    3 - fasta file for whole genome
    4 - bed file with bi-sulfide methylation direction and percentage
    
    Wren Saylor
    January 13 2017
    
    """

import argparse
import Bio
from Bio import SeqIO
from Bio import Seq
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from numpy import sin, linspace, pi
import numdifftools as nd
import numpy as np
import numpy.ma as ma
import pandas as pd
import pybedtools as pbt
import re
import scipy.fftpack
from scipy.interpolate import splrep, splev
from scipy import signal as ssignal
import scipy as sp
import scipy.stats as ss
from scipy.stats import mstats
import seaborn as sns
import tempfile

# set command line arguments
def get_args():
    parser = argparse.ArgumentParser(description="Description")
        parser.add_argument("efile", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files with unique names separated by newlines')
        # GENOME FILES
        parser.add_argument("mfile", type=argparse.FileType('rU'), help="A file containing a list of paths to the methylation files with unique names separated by newlines, with data for methylation position (chr, start,stop) and methylation % as fourth column'")
        parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
        #parser.add_argument("-n", "--nucleosome", type=str, help="A bedgraph file with data for nucleosome positions, form 'chr, start, stop, occupancy'")
        #parser.add_argument("-s", "--snp", type=str, help="A file with data for snps, form 'chr, start, stop(start+size alt-mutated), ref, ref_size, alt, alt_size, af_adj'")
        parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")
        return parser.parse_args()

# read in files
def eachFileProcess(fileName):
    btFeatures = pbt.BedTool(fileName)
        return btFeatures

# get bt features
def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
        return btFeatures

# get the correct range for fang evaluation
def getRange(btFeatures,fileName,num,uce,inuce):
    flankSize = (num - uce)/2
        inregion = uce-(inuce*2)
        midFeatures = pd.read_table(btFeatures.fn, header=None)
        midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
        midFeatures['sCenter'] = midFeatures['middle'] - (inregion/2)
        midFeatures['eCenter'] = midFeatures['middle'] + (inregion/2)
        midFeatures['sEdge'] = midFeatures.loc[:,1] + inuce
        midFeatures['eEdge'] = midFeatures.loc[:,2] - inuce
        midFeatures['sBoundary'] = midFeatures.loc[:,1] - flankSize
        midFeatures['eBoundary'] = midFeatures.loc[:,2] + flankSize
        midFeatures['start'] = midFeatures.loc[:,1]
        midFeatures['end'] = midFeatures.loc[:,2]
        midFeatures['type'] = midFeatures.loc[:,4]
        midFeatures['id'] = midFeatures.loc[:,3]
        midFeatures['chr'] = midFeatures.loc[:,0]
        midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
        rangeFeatures = midFeatures[['type','id','size','chr','sBoundary', 'start', 'sEdge', 'sCenter','eCenter','eEdge','end','eBoundary']] #,'starttwohund','endtwohund'
        return rangeFeatures

# save file from bedtool
def saveBedTool(btObject, strFilename):
    btObject.saveas(strFilename)

# convert bedtool to panda
# If there is nothing in the btobject, will it read the data from the previous itteration!?
def bedtoolToPanda(btobject):
    saveBedTool(btobject,'temp.bed')
        pdObject = pd.read_table(btobject.fn, header=None)
        return pdObject

# intersect the methylation file with the elements
def methIntersect(rangeFeatures,methFeature,num,uce,inuce):
    methSBoundary = methFeature.intersect(rangeFeatures[['chr','sBoundary','sEdge','id']].values.tolist(),wb=True,wa=True)
        pdmethSBoundary = bedtoolToPanda(methSBoundary)
        pdmethSBoundary['intAdd'] =  0
        pdmethSBoundary.columns = ['mchr','mstart','mstop','methPer','chr','sBoundary','sEdge','id','intAdd']
        pdmethSBoundary['methLoc'] = pdmethSBoundary['sEdge'].astype(int)-pdmethSBoundary['mstop'].astype(int)+pdmethSBoundary['intAdd'].astype(int) # or just count?
        outpdmethSBoundary = pdmethSBoundary[['id','methPer','methLoc']]
        methMiddle = methFeature.intersect(rangeFeatures[['chr','sCenter','eCenter','id']].values.tolist(),wb=True,wa=True)
        pdmethFeature = bedtoolToPanda(methMiddle)
        pdmethFeature['intAdd'] = (((num - uce)/2) + inuce)
        pdmethFeature.columns = ['mchr','mstart','mstop','methPer','chr','sCenter','eCenter','id','intAdd']
        pdmethFeature['methLoc'] = pdmethFeature['eCenter'].astype(int)-pdmethFeature['mstop'].astype(int)+pdmethFeature['intAdd'].astype(int) # or just count?
        outpdmethFeature = pdmethFeature[['id','methPer','methLoc']]
        methEBoundary = methFeature.intersect(rangeFeatures[['chr','eEdge','eBoundary','id']].values.tolist(),wb=True,wa=True)
        pdmethEBoundary = bedtoolToPanda(methEBoundary)
        pdmethEBoundary['intAdd'] = (((num - uce)/2) + (uce - inuce))
        pdmethEBoundary.columns = ['mchr','mstart','mstop','methPer','chr','eEdge','eBoundary','id','intAdd']#,'methDir'
        pdmethEBoundary['methLoc'] = pdmethEBoundary['eBoundary'].astype(int)-pdmethEBoundary['mstop'].astype(int)+pdmethEBoundary['intAdd'].astype(int) # or just count?
        outpdmethEBoundary = pdmethEBoundary[['id','methPer','methLoc']]
        methList = [outpdmethSBoundary,outpdmethFeature,outpdmethEBoundary]
        concatMeth = pd.concat(methList)
        #threshMeth = (concatMeth[(concatMeth['methPer'] >= 50)]) # subset for % methylation over 50
        sortMeth = concatMeth.sort_values(['methLoc'],ascending=True)
        #if want to collect by id, for use of cpg population http://stackoverflow.com/questions/10373660/converting-a-pandas-groupby-object-to-dataframe
        #groupMeth = sortMeth.groupby(['id'])['methLoc'].apply(list).to_frame(name = 'groupMeth').reset_index()
        return sortMeth

# collect methylation frequencies by position of elements
def methylationFreq(methPosition,num): # concatMeth from methIntersect
    new_index = range(0,num)
        subsetMeth = methPosition[['methLoc','methPer']]
        subsetMeth['methFreq'] = subsetMeth.groupby(['methLoc'])['methPer'].transform('count')
        dupMeth = subsetMeth.sort_values(['methLoc'],ascending=True).drop_duplicates(['methLoc','methFreq']) # returns highest value of methylation
        methIndex = dupMeth.set_index('methLoc').reindex(new_index,fill_value=0)
        methIndex.index.name = None
        intMeth = methIndex.astype(int)
        return intMeth

# get methylation positions for all methylation files
def methPositions(mFiles,rangeFeatures,num,uce,inuce):
    outMeth = []
        for methName in mFiles:
            methFeatures = eachFileProcess(methName)
                methPosition= methIntersect(rangeFeatures,methFeatures,num,uce,inuce)
                methFreq = methylationFreq(methPosition,num)
                methFreq.columns = ['{0}_Percentage'.format(methName),'{0}_Frequency'.format(methName)]
                outMeth.append(methFreq)
        pdMeth = pd.concat(outMeth,axis=1)
    return pdMeth

# get fasta strings for each desired region
def getFasta(btFeatures,faGenome,fileName):
    saveAll = 'Seq_result_for_all_{0}.txt'.format(fileName)
        seqAll = btFeatures.sequence(fi=faGenome,fo=saveAll)
        saveExonic = 'Seq_results_for_exonic_{0}.txt'.format(fileName)
        seqExonic = btFeatures.sequence(fi=faGenome,fo=saveExonic).filter(lambda x: x[name] == 'exonic')
        saveIntronic = 'Seq_results_for_intronic_{0}.txt'.format(fileName)
        seqIntronic = btFeatures.sequence(fi=faGenome,fo=saveIntronic).filter(lambda x: x[name] == 'intronic')
        saveIntergenic = 'Seq_results_for_intergenic_{0}.txt'.format(fileName)
        seqIntergenic = btFeatures.sequence(fi=faGenome,fo=saveIntergenic).filter(lambda x: x[name] == 'intergenic')
        return saveAll, saveExonic, saveIntronic, saveIntergenic

# used in btRange to extract just the fasta strings
def simpleFasta(inFeature,faGenome):
    seqFeature = inFeature.sequence(fi=faGenome)
        outFeature = pd.read_table(seqFeature.seqfn)
        outSequence = outFeature[::2]
        outSequence = outSequence.reset_index(drop=True)
        return outSequence

# get the strings for sliding window regions
def btRange(rangeFeatures,faGenome):#,methFeature
    rangeFeatures['sBoundarySeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','sBoundary','start']].values.tolist()),faGenome)
        rangeFeatures['sEdgeSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','start','sEdge']].values.tolist()),faGenome)
        rangeFeatures['MiddleSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','sCenter','eCenter']].values.tolist()),faGenome)
        rangeFeatures['eEdgeSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','eEdge','end',]].values.tolist()),faGenome)
        rangeFeatures['eBoundarySeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','end','eBoundary']].values.tolist()),faGenome)
        rangeFeatures['feature'] = simpleFasta(getFeatures(rangeFeatures[['chr','start','end']].values.tolist()),faGenome)
        rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
        rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
        rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
        rangeFeatures['reverseComplement'] = rangeFeatures.apply(lambda row: reverseComplement(row['combineString']),axis=1)
        # 	methPosition = methIntersect(rangeFeatures,methFeature)
        # 	outFeatures = pd.merge(rangeFeatures,methPosition,left_on='id',right_on='id',how='outer',indicator=True) #perhaps want to change indicator and how args
        # 	outFeatures['groupMeth'].fillna('[-999]',inplace=True) #,-999 have to trick into accepting a list, for iterating through later
        return rangeFeatures

# sliding window
def slidingWindow(element,num,uce,inuce,window):
    winFeatures = []
        winCpG = []
        winA = []
        winT = []
        winG = []
        winC = []
        winMotif = []
        n = num #600 # len(element) # take out hard coding
        start, end = 0, window
        while end < n:
            current = element[start:end]
                #print num, float(len(element)), float(len(current)), start, end, current
                percentageAT = eval('100 * float(current.count("A") + current.count("T"))/ float(len(current))')
                percentageCpG = eval('100 * float(current.count("CG")) / float(len(current))')
                perA = eval('100 * float(current.count("A"))/ float(len(current))')
                perT = eval('100 * float(current.count("T"))/ float(len(current))')
                perG = eval('100 * float(current.count("G"))/ float(len(current))')
                perC = eval('100 * float(current.count("C"))/ float(len(current))')
                perMo = eval('100 * float(current.count("ATTAAT")) / float(len(current))')
                winFeatures.append(percentageAT)
                winCpG.append(percentageCpG)
                winA.append(perA)
                winT.append(perT)
                winG.append(perG)
                winC.append(perC)
                winMotif.append(perMo)
                start, end = start + 1, end + 1
        return winFeatures, winCpG, winA, winT, winG, winC, winMotif

# append results from sliding window
def appendWindow(rangeFeatures,num,uce,inuce,window):
    outWindow = []
        outCpG = []
        outA = []
        outT = []
        outG = []
        outC = []
        outMo = []
        for element in rangeFeatures['combineString']:
            outFeature, winCpG, winA, winT, winG, winC, winMo = slidingWindow(element,num,uce,inuce,window)
                outWindow.append(outFeature)
                outCpG.append(winCpG)
                outA.append(winA)
                outT.append(winT)
                outG.append(winG)
                outC.append(winC)
                outMo.append(winMo)
        return outWindow, outCpG, outA, outT, outG, outC, outMo

# run append to sliding window and return pandas data frames
def dataframeWindow(rangeFeatures,num,uce,inuce,window):
    outWindow, outCpG, outA, outT, outG, outC, outMo = appendWindow(rangeFeatures,num,uce,inuce,window)
        pdWindow = pd.DataFrame(outWindow)
        pdCpG = pd.DataFrame(outCpG)
        pdA = pd.DataFrame(outA)
        pdT = pd.DataFrame(outT)
        pdG = pd.DataFrame(outG)
        pdC = pd.DataFrame(outC)
        pdMo = pd.DataFrame(outMo)
        return pdWindow, pdCpG, pdA, pdT, pdG, pdC, pdMo

# Collect information about how each UCEs second derivative behaves
def behaviorUCE(fillX,pdWindow):
    for index, row in pdWindow.iterrows():
        f = splrep(fillX,row,k=5,s=11)
            smoothMean = splev(fillX,f)
                secondDer = splev(fillX,f,der=2)
                secondDer[0:window] = 0 # small edge effect
                secondDer[-window:] = 0 # small edge effect
                return secondDer

# make some graphs!
def endLinegraphs(pdMeth,pdWindow,pdCpG, pdA,pdT,pdG,pdC,pdMo,fileName,num,uce,inuce,window):
    fillX = range(0,(num-window))
        halfwindow = ((window/2)+1)
        sns.set_style('ticks')
        gs = gridspec.GridSpec(3,3,height_ratios=[3,1,1]) # for the ratios of the graphs against each other
        gs.update(hspace=.8) # setting the space between the graphs
        info = str(fileName) + ', '+ str(len(pdWindow)) + ' - ' "UCES"
        plt.suptitle(info,fontsize=10)
        mean = pdWindow.mean()
        #maX = np.ma.masked_where(240<fillX<250 & 340<fillX<350, fillX) # tried to make a mast for gapped data
        pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
        
        # Plot the mean AT content with a std of 1
        # NEED TO USE THE RATE CHANGE GRAPH TO GET THE LOCATIONS TO FIND FANG!! .diff()
        # 	maxStartY = int(mean[150:250].max())
        # 	maxStopY = int(mean[350:450].max())
        # 	minStartY = int(mean[150:250].min())
        # 	minStopY = int(mean[350:450].min())
        # 	maxStartX = fillX[mean[150:250].argmax(maxStartY)]
        # 	maxStopX = fillX[mean[350:450].argmax(maxStopY)]
        # 	minStartX = fillX[mean[150:250].argmin(minStartY)]
        # 	minStopX = fillX[mean[350:450].argmin(minStopY)]
        StartMean = pdWindow.loc[:,(((num-uce)/2)-halfwindow-window):(((num-uce)/2)+(inuce-halfwindow))].mean(axis=0) # 185*:245 - trying to just have the boundary and 50inward (184)
        StopMean = pdWindow.loc[:,(((num-uce)/2)+(uce-inuce-halfwindow)):(((num-uce)/2)+uce-(halfwindow-window))].mean(axis=0) # 345:*405 - trying to just have the boundary and 50inward (406)
        #StopMean.iloc[::-1]
        wilcoxPSRMean = ss.wilcoxon(StartMean,StopMean)
        ax0 = plt.subplot(gs[0,:])
        ax0.plot(fillX,mean,linewidth=1, color='#3e1638',label='AT')
        ax0.fill_between(fillX,mean+pdWindow.std(axis=0,ddof=1),mean-pdWindow.std(axis=0,ddof=1),facecolor = '#63245a',label='',alpha=0.3)
        ax0.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a') # 245
        ax0.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a') # 345
        ax0.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973') # 195
        ax0.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973') # 395
        ax0.hlines(y=86,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
        ax0.text(32,85,'11bp sliding window',size=6)
        ax0.text(20,90,'Wilcox Signed Rank P-value {:0.1e}'.format(wilcoxPSRMean[1]),size=6,clip_on=False)
        ax0.set_ylabel('% AT Content',size=8)
        ax0.set_xlabel('Position',size=6)
        ax0.legend(loc=0,fontsize=5,labelspacing=0.1)
        ax0.set_title('Mean AT Content With Standard Deviation',size=8)
        ax0.set_yticks(ax0.get_yticks()[::2])
        plt.xlim(0,num)
        
        # Plot the std = 1
        ax1 = plt.subplot(gs[1,:],sharex=ax0)
        ax1.plot(fillX,pdWindow.std(axis=0,ddof=1),linewidth=1, color='#3e1638',label='AT')
        ax1.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax1.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax1.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax1.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax1.axvspan(0,(((num-uce)/2)-halfwindow),facecolor = '#863eae',label='',alpha=0.1) # 0,195
        ax1.axvspan((((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow),facecolor = '#ae3e9e',label='',alpha=0.1) #195,395
        ax1.axvspan((((num-uce)/2)+uce-halfwindow),(num-window),facecolor = '#ae3e66',label='',alpha=0.1) # 395,589
        ax1.set_yticks(ax1.get_yticks()[::2])
        ax1.set_xlabel('Position',size=6)
        ax1.set_ylabel('SD',size=8)
        ax1.set_title('Standard Deviation',size=8) #, With One Degree of Freedom
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.legend(loc=0,fontsize=5,labelspacing=0.05)
        
        # Significances tests for SD populations
        upStream = pdWindow.loc[:,0:(((num-uce)/2)-halfwindow)].std(axis=0,ddof=1) # 0,195
        downStream = pdWindow.loc[:,(((num-uce)/2)+uce-halfwindow):(num-window)].std(axis=0,ddof=1) # 395,589
        uceRegion = pdWindow.loc[:,(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)].std(axis=0,ddof=1) ##195,395
        bothStream = pdWindow.iloc[:,np.r_[0:(((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow):(num-window)]].std(axis=0,ddof=1) # 0,195,395,589
        # Kruskal-Wallis test
        # 	kruskalUP = mstats.kruskalwallis(upStream,uceRegion)#$p.value
        # 	kruskalDown = mstats.kruskalwallis(downStream,uceRegion)
        # 	kruskalBoth = mstats.kruskalwallis(upStream,downStream)
        # 	dfKruskal = pd.DataFrame(['{:0.2e}'.format(kruskalBoth[1]),'{:0.2e}'.format(kruskalUP[1]),'{:0.2e}'.format(kruskalDown[1])],index =['Up - Down','Up - UCE','Down - UCE'])
        kruskalSD = mstats.kruskalwallis(bothStream,uceRegion)
        ax2 = plt.subplot(gs[2,0])
        ax2.hist(upStream,35,linewidth=0.3, color='#3e1638',label='UpStream',alpha=0.7)
        ax2.set_yticks(ax2.get_yticks()[::2])
        ax2.set_ylabel('Frequency',size=8)
        ax2.set_xlabel('SD Value',size=6)
        ax3 = plt.subplot(gs[2,1],sharey=ax2)
        # 	ax3.text(16.75,11,'B P-value {:0.1e}'.format(bartlettSD[1]),size=6)
        # 	ax5.text(16.75,17,'KS P-value {:0.1e}'.format(kolmogorovSD[1]),size=6,clip_on=False)
        # 	ax5.text(16.75,14,'AD P-value {:0.1e}'.format(andersonSD[2]),size=6,clip_on=False)
        ax3.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        ax3.hist(uceRegion,35,linewidth=0.3, color='#ae3e9e',label='UpStream',alpha=0.7)
        ax3.set_title('Standard Deviation Frequency for Highlighted Regions',size=8)
        ax3.set_xlabel('SD Value',size=6)
        ax4 = plt.subplot(gs[2,2],sharey=ax2)
        plt.setp(ax4.get_yticklabels(), visible=False)
        ax4.hist(downStream,35,linewidth=0.3, color='#ae3e66',label='UpStream',alpha=0.7)
        ax4.set_xlabel('SD Value',size=6)
        # Table for Kruskal Test P-value results
        # 	ax7 = plt.subplot(gs[3,3])
        # 	ax7.set_frame_on(False)
        # 	ax7.set_yticks([])
        # 	ax7.set_xticks([])
        # 	KWTable = ax7.table(cellText=dfKruskal.values,rowLabels=dfKruskal.index,cellLoc='center',rowLoc='center',loc='center right',colWidths=(.5,.5))
        # 	KWTable.set_fontsize(6)
        # 	ax7.set_title('P-values for Kruskal-Wallis',size=8)
        sns.despine()
        plt.savefig(pp, format='pdf')
        
        gs = gridspec.GridSpec(4,3,height_ratios=[1,1,1,1])
        gs.update(hspace=.65)
        
        # Plot mean AT
        ax5 = plt.subplot(gs[0,:],sharex=ax0)
        ax5.plot(fillX,pdWindow.mean(axis=0),linewidth=1, color='#3e1638')
        ax5.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax5.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax5.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax5.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax5.hlines(y=63,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
        ax5.text(32,63,'11bp sliding window',size=6)
        ax5.set_yticks(ax5.get_yticks()[::2])
        ax5.set_ylabel('% AT Content',size=8)
        ax5.set_title('Mean AT Content',size=8)
        
        # Create fitted, first and second derivative lines
        f = splrep(fillX,pdWindow.mean(axis=0),k=5,s=11)
        smoothMean = splev(fillX,f)
        firstDer = splev(fillX,f,der=1)
        firstDer[0:halfwindow] = 0 # small edge effect
        firstDer[-halfwindow:] = 0 # small edge effect
        secondDer = splev(fillX,f,der=2)
        secondDer[0:window] = 0 # small edge effect
        secondDer[-window:] = 0 # small edge effect
        
        # Plot smoothed mean AT
        ax6 = plt.subplot(gs[1,:],sharex=ax0)
        ax6.plot(fillX,smoothMean,linewidth=1, color='#3e1638',alpha=0.9)
        ax6.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax6.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax6.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax6.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax6.set_yticks(ax6.get_yticks()[::2])
        ax6.set_ylabel('% AT Content',size=8)
        ax6.set_title('Fitted Mean AT Content',size=8)
        
        # First derivative
        ax7 = plt.subplot(gs[2,:],sharex=ax0)
        ax7.plot(fillX,firstDer,linewidth=1, color='#3e1638',alpha=0.8) # mean.diff(), np.gradient(mean), numpy.diff(y) / numpy.diff(x)
        ax7.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax7.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax7.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax7.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax7.axhline(y=0,linewidth=.1,color='#bd4973',alpha=0.3)
        ax7.set_yticks(ax7.get_yticks()[::2])
        ax7.set_ylabel('Amplitude',size=8)
        ax7.set_title('First Derivative of Fitted Mean',size=8)
        
        # Second derivative
        ax8 = plt.subplot(gs[3,:],sharex=ax0)
        ax8.plot(fillX,secondDer,linewidth=1, color='#3e1638',alpha=0.7)
        #http://stackoverflow.com/questions/13691775/python-pinpointing-the-linear-part-of-a-slope
        #http://stackoverflow.com/questions/16323139/finding-inflection-points-in-spline-fitted-1d-data
        ax8.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax8.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax8.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax8.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax8.axhline(y=0,linewidth=.1,color='#bd4973',alpha=0.3)
        ax8.set_ylabel('Amplitude',size=8)
        ax8.set_xlabel('Position',size=6)
        ax8.set_yticks(ax8.get_yticks()[::2])
        ax8.set_title('Second Derivative of Fitted Mean',size=8)
        sns.despine()
        pp.savefig()
        
        gs = gridspec.GridSpec(3,3,height_ratios=[2,1,1])
        gs.update(hspace=.65)
        
        # Short Fourier Transform
        ax9 = plt.subplot(gs[0,:],sharex=ax0)
        f1, t1, Zxx1 = ssignal.stft(firstDer,fs=1.0, window='hann',nperseg=30,noverlap=None)#,nperseg=11,noverlap=5
        ax9.pcolormesh(t1,f1,np.abs(Zxx1),cmap='RdPu')
        ax9.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a') # 245
        ax9.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a') # 345
        ax9.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973') # 195
        ax9.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973') # 395
        ax9.set_ylabel('Frequency',size=8)
        ax9.set_xlabel('Position',size=6)
        ax9.set_yticks(ax9.get_yticks()[::2])
        ax9.set_title('Short Fourier Transform',size=8)#30 bp bins
        
        # First Derivative
        ax10 = plt.subplot(gs[1,:],sharex=ax0)
        ax10.plot(fillX,firstDer,linewidth=1, color='#3e1638')
        ax10.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax10.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax10.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax10.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax10.axvspan((((num-uce)/2)+(inuce-halfwindow)),(((num-uce)/2)+(uce-inuce-halfwindow)),facecolor = '#ae3e9e',label='',alpha=0.1) # 245, 345
        ax10.axvspan(inuce,(((num-uce)/2)-inuce),facecolor = '#863eae',label='',alpha=0.1) # 50, 150
        ax10.axvspan((num-window-(((num-uce)/2)-inuce)),(num-window-inuce),facecolor = '#ae3e66',label='',alpha=0.1) # 439,539
        ax10.set_yticks(ax10.get_yticks()[::2])
        ax10.set_xlabel('Position',size=6)
        ax10.set_ylabel('Amplitude',size=8)
        ax10.set_title('First Derivative of Fitted Mean',size=8)
        
        #https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.signal.stft.html
        #https://mathematica.stackexchange.com/questions/18082/plotting-the-frequency-spectrum-of-a-data-series-using-fourier
        #http://stackoverflow.com/questions/1523814/units-of-a-fourier-transform-fft-when-doing-spectral-analysis-of-a-signal
        Fs = 1.0 # sampling rate
        Ts = 1.0/Fs # sampling interval
        y2sd = firstDer[(((num-uce)/2)+(inuce-halfwindow)):(((num-uce)/2)+(uce-inuce-halfwindow))] # 245, 345
        n2sd = len(y2sd) # length of the signal
        k2sd = np.arange(n2sd)
        T2sd = n2sd/Fs
        frq2sd = k2sd/T2sd # two sides frequency range
        frq2sd = frq2sd[range(n2sd/2)] # one side frequncy range
        Y2sd = np.fft.fft(y2sd)/n2sd # fft computing and normalization
        Y2sd = Y2sd[range(n2sd/2)]
        y3sd = firstDer[inuce:(((num-uce)/2)-inuce)] # 50, 150
        n3sd = len(y3sd)
        k3sd = np.arange(n3sd)
        T3sd = n3sd/Fs
        frq3sd = k3sd/T3sd
        frq3sd = frq3sd[range(n3sd/2)]
        Y3sd = np.fft.fft(y3sd)/n3sd
        Y3sd = Y3sd[range(n3sd/2)]
        y4sd = firstDer[(num-window-(((num-uce)/2)-inuce)):(num-window-inuce)]#439,539
        n4sd = len(y4sd)
        k4sd = np.arange(n4sd)
        T4sd = n4sd/Fs
        frq4sd = k4sd/T4sd
        frq4sd = frq4sd[range(n4sd/2)]
        Y4sd = np.fft.fft(y4sd)/n4sd
        Y4sd = Y4sd[range(n4sd/2)]
        
        # 	FFT for sections of the smoothed second derivative
        ax11 = plt.subplot(gs[2,0])
        ax11.plot(frq3sd,abs(Y3sd),linewidth=1, color='#863eae')
        ax11.set_ylabel('|Y(freq)|',size=8)
        ax11.set_xlabel('Freq(Hz)',size=6)#AT Rate Change
        ax11.set_yticks(ax11.get_yticks()[::2])
        ax12 = plt.subplot(gs[2,1],sharey=ax11)
        plt.setp(ax12.get_yticklabels(), visible=False)
        ax12.plot(frq2sd,abs(Y2sd),linewidth=1, color='#ae3e9e')
        ax12.set_title('Power Series for Highlighted Regions',size=8)# Power Spectrum Analysis for FFT
        ax12.set_xlabel('Freq(Hz)',size=6)
        ax13 = plt.subplot(gs[2,2],sharey=ax11)
        plt.setp(ax13.get_yticklabels(), visible=False)
        ax13.plot(frq4sd,abs(Y4sd),linewidth=1, color='#ae3e66')
        ax13.set_xlabel('Freq(Hz)',size=6)
        sns.despine()
        pp.savefig()
        
        gs = gridspec.GridSpec(2,1,height_ratios=[1,3])
        gs.update(hspace=.45)
        
        # Nucleotides
        ax14 = plt.subplot(gs[0],sharex=ax0)
        ax14.plot(fillX,pdA.mean(axis=0),linewidth=1, color='#3f1bd7',label='A')
        ax14.plot(fillX,pdT.mean(axis=0),linewidth=1, color='#d7401b',label='T')
        ax14.plot(fillX,pdG.mean(axis=0),linewidth=1, color='#d73f1b',label='G')
        ax14.plot(fillX,pdC.mean(axis=0),linewidth=1, color='#9d1bd7',label='C')
        ax14.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax14.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax14.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax14.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax14.set_yticks(ax14.get_yticks()[::2])
        ax14.set_ylabel('% Nucleotide Content',size=8)
        ax14.set_xlabel('Position',size=6)
        ax14.set_title('Mean Nucleotide Content Separately',size=8)
        ax14.legend(loc=0,fontsize=5,labelspacing=0.05)
        
        # Plot the CpG
        # Might still want to return the actual CpG location for how many are methylated
        ax15 = plt.subplot(gs[1],sharex=ax0)
        meanCpG = pdCpG.mean(axis=0)
        ax15.plot(fillX,meanCpG,linewidth=1, color='#d71b54',label='CpG mean',alpha=0.7) # orange
        ax15.fill_between(fillX,pdCpG.mean(axis=0)+pdCpG.std(axis=0,ddof=1),pdCpG.mean(axis=0)-pdCpG.std(axis=0,ddof=1),facecolor = '#d71b54',label='',alpha=0.2)
        ax15.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax15.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax15.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax15.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax15.set_title('Mean CpG',size=8)
        ax15.set_ylabel('% CpG Content',size=8)
        ax15.set_xlabel('Position',size=6)
        ax15.set_yticks(ax15.get_yticks()[::2])
        sns.despine()
        pp.savefig()
        
        gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
        gs.update(hspace=.5)
        pdMethPer = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('Percentage',case=False)]])
        pdMethNum = (pdMeth[pdMeth.columns[pdMeth.columns.str.contains('Frequency',case=False)]])
        TPer = pdMethPer.T
        TNum = pdMethNum.T
        
        # Make heatmap for % methylation (Percentage)
        ax16 = plt.subplot(gs[0,:],sharex=ax0)
        heatmap1 = ax16.pcolormesh(TPer,cmap='RdPu')#col_cluster=False sns.clustermap
        ax16.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax16.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax16.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax16.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax16.set_ylabel('% Methylation per Sample',size=8)
        ax16.set_xlabel('Position',size=6)
        ax16.tick_params(labelsize=8)
        ylabels1 = TPer.index.str.replace('.bed_Percentage','')
        ax16.set_yticklabels(ylabels1,minor=False)
        ax16.set_yticks(np.arange(TPer.shape[0]) + 0.5, minor=False)
        ax16.set_title('Methylation Percentage',size=8)
        cbar1 = plt.colorbar(mappable=heatmap1,orientation="vertical",shrink=.7,pad=0.072)#,label="% Methylation" ,anchor=(0.0, 0.5)
        cbar1.set_clim(vmin=0, vmax=100)
        
        # Make heatmap for # methylation (Frequency)
        ax17 = plt.subplot(gs[1,:],sharex=ax0)
        heatmap2 = ax17.pcolormesh(TNum,cmap='RdPu')#col_cluster=False sns.clustermap
        ax17.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax17.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#e7298a')
        ax17.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax17.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#bd4973')
        ax17.set_ylabel('Frequency of Methylation per Sample',size=8)
        ax17.set_xlabel('Position',size=6)
        ax17.tick_params(labelsize=8)
        ylabels2 = TNum.index.str.replace('.bed_Frequency','')
        ax17.set_yticklabels(ylabels2,minor=False)
        ax17.set_yticks(np.arange(TNum.shape[0]) + 0.5, minor=False)
        ax17.set_title('Methylation Frequency',size=8)
        cbar2 = plt.colorbar(mappable=heatmap2,orientation="vertical",shrink=.7,pad=0.072)#,label="% Methylation" ,anchor=(0.0, 0.5)
        cbar2.set_clim(vmin=0, vmax=5)
        
        # Heat map for interaction of % and #
        
        # Heat map for % methlation of cpg
        
        sns.despine()
        pp.savefig()
        pp.close()

# get the reverse complement
def reverseComplement(sequence):
    seqDict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join([seqDict[base] for base in reversed(sequence)])

# separate by directionality
def dirLine(rangeFeatures,fileName,mFiles,num,uce,inuce,window):
    DirList = ["+","-"]#,"="
        for direction in DirList:
            dirStr = (rangeFeatures[(rangeFeatures['compareBoundaries'] == direction)])
                #savePanda(dirStr[['chr','start','end','type']], '{0}_30_coordinates_{1}.bed'.format(direction,fileName)) # prints the bed file for each directionality
                dirWindow, dirWinCpG, dirWinA, dirWinT, dirWinG, dirWinC, dirWinMo = dataframeWindow(dirStr,num,uce,inuce,window)
                pddirMeth = methPositions(mFiles,dirStr,num,uce,inuce)
                endLinegraphs(pddirMeth,dirWindow,dirWinCpG,dirWinA,dirWinT,dirWinG,dirWinC,dirWinMo,'{0}_30_{1}'.format(direction,fileName),num,uce,inuce,window)
        # More drawn out performance of sliding window and methylation collection in order to collect reverse complement super imposed
        negStr = (rangeFeatures[(rangeFeatures['compareBoundaries'] == '-')])
        outComp = []
        outCpG = []
        compA = []
        compT = []
        compG = []
        compC = []
        compMo = []
        outnegMeth = []
        for element in negStr['reverseComplement']:
            negFeature, negCpGfeature,negAfeature,negTfeature,negGfeature,negCfeature,negMofeature = slidingWindow(element,num,uce,inuce,window)
                outComp.append(negFeature)
                outCpG.append(negCpGfeature)
                compA.append(negAfeature)
                compT.append(negTfeature)
                compG.append(negGfeature)
                compC.append(negCfeature)
                compMo.append(negMofeature)
        for methName in mFiles:
            methFeatures = eachFileProcess(methName)
                methPosition = methIntersect(negStr,methFeatures,num,uce,inuce)
                methFreq = methylationFreq(methPosition,num)
                methFreq.columns = ['{0}_Percentage'.format(methName),'{0}_Frequency'.format(methName)]
                outnegMeth.append(methFreq)
posStr = (rangeFeatures[(rangeFeatures['compareBoundaries'] == '+')])
    for element in posStr['combineString']:
        posFeature, posCpGfeature,posAfeature, posTfeature, posGfeature, posCfeature,posMofeature = slidingWindow(element,num,uce,inuce,window)
            outComp.append(posFeature)
                outCpG.append(posCpGfeature)
                compA.append(posAfeature)
                compT.append(posTfeature)
                compG.append(posGfeature)
                compC.append(posCfeature)
                compMo.append(posMofeature)
        compWindow = pd.DataFrame(outComp)
        compWinCpG = pd.DataFrame(outCpG)
        compWinA = pd.DataFrame(compA)
        compWinT = pd.DataFrame(compT)
        compWinG = pd.DataFrame(compG)
        compWinC = pd.DataFrame(compC)
        compWinMo = pd.DataFrame(compMo)
        pdposMeth = methPositions(mFiles,posStr,num,uce,inuce)
        pdnegMeth = pd.concat(outnegMeth,axis=1)
        pdnegInMeth = pdnegMeth.iloc[::-1] # reverse order 
        pdnegInsetMeth = pdnegInMeth.reset_index(drop=True) # make new index, don't make old one into column
        pdcompMeth = pd.concat([pdnegInsetMeth,pdposMeth],axis=1)
        pdcompMeth = pdcompMeth.groupby(pdcompMeth.columns, axis=1).sum()# sum-Freq, max-Perc
        # 	pdcompMeth = pdcompMeth.groupby(pdcompMeth.columns, axis=1).agg({pdcompMeth.columns.str.contains('Frequency',case=False):np.sum,pdcompMeth.columns.str.contains('Percentage',case=False):lambda x: np.max(x)})
        # 	print pdcompMeth
    endLinegraphs(pdcompMeth,compWindow,compWinCpG,compWinA,compWinT,compWinG,compWinC,compWinMo,'revComp_30_{0}'.format(fileName),num,uce,inuce,window)

# save file from panda
def savePanda(pdData, strFilename):
    pdData.to_csv(strFilename, sep='\t', header=False, index=False)

# out put directionality, as inferred by comparing first and last n base pairs
def compareN(element,size):
    start = element[:size]
        end = element[-size:]
        perSize = []
        perSize.append(eval('100*int(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
        perSize.append(eval('100*int(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
        if perSize[0] > perSize[1]: outList = '+'
        if perSize[1] > perSize[0]: outList = '-'
        if perSize[1] == perSize[0]: outList = '='
        return outList

# with the results from compareN per each element, evaluate directionality into new column
def evalN(rangeFeatures,fileName,binDir):
    rangeFeatures['compareBoundaries'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],binDir)),axis=1)
        compareEnds = pd.DataFrame(rangeFeatures[['chr','start','end','compareBoundaries']])
        return rangeFeatures

# do all the analysis/plots for each type
def perType(rangeFeatures,fileName,mFiles,num,uce,inuce,window):
    typeList = ['intergenic','intronic','exonic']
        for type in typeList:
            boolType = (rangeFeatures[rangeFeatures['type'] == type])
                pdWindow, pdTypeCpG, pdTypeA, pdTypeT, pdTypeG, pdTypeC, pdTypeMo = dataframeWindow(boolType,num,uce,inuce,window)
                pdMeth= methPositions(mFiles,boolType,num,uce,inuce)
                endLinegraphs(pdMeth,pdWindow,pdTypeCpG,pdTypeA,pdTypeT,pdTypeG,pdTypeC,pdTypeMo,'{0}_{1}'.format(type,fileName),num,uce,inuce,window)
                dirLine(boolType,'{0}_{1}'.format(type,fileName),mFiles,num,uce,inuce,window)

def main():
    # Collect arguments
    num = 600 # 600 # total size of region to look at (region + flanks), even number # suggested to be at least triple your element
        uce = 200 # size of your element (region without flanks), even number
        inuce = 50 # size into your element from the boundaries, even number, suggested 50
        window = 11 # size of sliding window, odd number, suggested 11
        binDir = 30 # size of bins to compare element ends, suggested 30
        args = get_args()
        eFiles = [line.strip() for line in args.efile]
        mFiles = [line.strip() for line in args.mfile]
        sizeGenome = args.genome
        faGenome = args.fasta
        #snpFeature = eachFileProcess(args.snp)
        for fileName in eFiles:
            btFeatures = eachFileProcess(fileName)
                subsetFeatures = getRange(btFeatures, fileName,num,uce,inuce)
                #saveAll, saveExonic, saveIntronic, saveIntergenic = getFasta(btFeatures,faGenome,fileName)
                rangeFeatures = btRange(subsetFeatures,faGenome)
                #http://stackoverflow.com/questions/19828822/how-to-check-whether-a-pandas-dataframe-is-empty
                pdWindow, pdCpG, pdA, pdT, pdG, pdC, pdMo = dataframeWindow(rangeFeatures,num,uce,inuce,window)
                pdMeth = methPositions(mFiles,rangeFeatures,num,uce,inuce)
                endLinegraphs(pdMeth,pdWindow,pdCpG,pdA,pdT,pdG,pdC,pdMo,fileName,num,uce,inuce,window)
# 		directionFeatures = evalN(rangeFeatures,fileName,binDir)
# 		dirLine(directionFeatures,fileName,mFiles,num,uce,inuce,window)
# 		perType(directionFeatures,fileName,mFiles,num,uce,inuce,window)

if __name__ == "__main__":
    main()
