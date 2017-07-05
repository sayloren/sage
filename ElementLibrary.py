"""
Script to process the elements you what to analyze

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
import tempfile

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
	return rangeFeatures

# get the reverse complement
def reverseComplement(sequence):
	seqDict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# save file from bedtool
def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

# convert bedtool to panda
# If there is nothing in the btobject, will it read the data from the previous itteration!?
def bedtoolToPanda(btobject):
	saveBedTool(btobject,'temp.bed')
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

# convert panda to bedtool
def pandaToBedtool(panda): # columns 5+
	arArFeatures = panda.values.tolist()
	btoutFeatures = getFeatures(arArFeatures)
	return btoutFeatures

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

def main(num,uce,inuce,window,binDir,fileName,sizeGenome,faGenome):
	btFeatures = eachFileProcess(fileName)
	subsetFeatures = getRange(btFeatures, fileName,num,uce,inuce)
	rangeFeatures = btRange(subsetFeatures,faGenome)
	return rangeFeatures

if __name__ == "__main__":
	main()