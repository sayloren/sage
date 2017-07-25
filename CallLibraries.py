"""
Script to call and execute different sets of analyses for Element Structure

Wren Saylor
July 5 2017

To Do:
IO - Have input methylation data in separate callable folder, print output graphs to new folder
If no methylation plotting, don't process methylation files
Add verbose as a setting to args
Cluster plot - just the element for clustering, return out clustered info for organization, label row colors and column colors with tissue for id x pos
Graph for bin = size calculation
Table for each group % =
Random regions printed on graph
AT mean cluster SD row color

"""

import argparse
import ElementLibrary
import DirectionLibrary
import FangsLibrary
import MethylationLibrary
import RevCompLibrary
import TypeLibrary
import BinLibrary
import GraphFangLibrary
import GraphMethLibrary
import GraphSignalLibrary
import BokehLibrary
import GraphTableLibrary
import GraphClusterLibrary
import os

# set command line arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
# 	parser.add_argument('-d',help='directory to look for external files')
	parser.add_argument("efile", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("mfile", type=argparse.FileType('rU'), help="A file containing a list of paths to the methylation files with unique names separated by newlines, with data for methylation position (chr, start,stop) and methylation % as fourth column'")
# 	parser.add_argument("rfile",type=argparse.FileType('rU'), help="A file containing a list of paths to the random regions equable with your elements to plot in contrast")

	# Genome Files
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
# 	parser.add_argument("-n", "--nucleosome", type=str, help="A bedgraph file with data for nucleosome positions, form 'chr, start, stop, occupancy'")
# 	parser.add_argument("-s", "--snp", type=str, help="A file with data for snps, form 'chr, start, stop(start+size alt-mutated), ref, ref_size, alt, alt_size, af_adj'")
	parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")

	# Integer Parameters
	parser.add_argument("-t", "--total", type=int, default="600", help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e", "--element", type=int, default="200", help='size of your element (region without flanks), should be an even number')
	parser.add_argument("-i", "--inset", type=int, default="50", help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w", "--window", type=int, default="11", help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b", "--bin", type=int, default="30", help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-mc", "--thresholdcoverage", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')
	parser.add_argument("-mp", "--thresholdpercentage", type=int, default="0", help='size to threshold % methylation data')
	parser.add_argument("-mf", "--methylationflank", type=int, default="20", help='The number of base pairs to look at outside of the element for the methylation clusterplots')

	# Directionality Calibration
	parser.add_argument("-n", "--base", type=int, default=2)
	parser.add_argument("-c", "--combinations", type=int, default=30)

	# Specify which groups and graphs to run
	parser.add_argument('-type', "--elementype", default=[], nargs='*', choices=['all','intronic','exonic','intergenic'],help='which group types of element to run') # type
	parser.add_argument('-dir', "--elementdirection", default=[], nargs='*', choices=['+','-','='], help='which group direction of element to run') # direction
	parser.add_argument('-rc', "--reversecomplement",action='store_true', help='if reverse complement sorting required')
	parser.add_argument('-p',"--plots",default=[],nargs='*',choices=['fangs','methylation','signal','interactive','cluster'],help='the available graphs to plot') # 'table', 'combinations'
	parser.add_argument('-nuc',"--nucleotideline",default=['A','T'],nargs='+',help='type the nucleotide string combinations to search for in the element')#'CA','CT','CC','CG'
	parser.add_argument('-str',"--stringname",type=str,help='string to add to the outfile name')

	# Add additional descriptive file name information
# 	parser.add_argument('-v',help='print statements',store=True)
	return parser.parse_args()

def getFilepaths(path):
	relative_path = path
	current_dir = os.getcwd()
	return os.join(current_dir,relative_path)

# for type and direction, separate the groups and run the analyses
def groupSeparate(List,directionFeatures,typecolumn,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):

	# subset by bool presence
	bool = (directionFeatures[directionFeatures[typecolumn] == List])

	# if there is nothing in that set, skip
	if len(bool.index) != 0:
		Meth,Window,Names = TypeLibrary.main(bool,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
	return bool,Meth,Window,Names

# the plotting options, if in the list of plot flags, run graph
def plotGraphs(pdMeth,slidingWinDF,names,fileName,num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank):
	if 'fangs' in graphs:
		GraphFangLibrary.main(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine)
	if 'signal' in graphs:
		GraphSignalLibrary.main(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine)
	if 'methylation' in graphs:
		GraphMethLibrary.main(slidingWinDF,pdMeth,fileName,num,uce,inuce,window)
	if 'interactive' in graphs:
		BokehLibrary.main(slidingWinDF,fileName,num,uce,inuce,window,nucLine)
# 	if 'table' in graphs:
# 		GraphTableLibrary.main(slidingWinDF,fileName,num,uce,inuce,window)
# 	if 'combinations' in graphs:
# 		BinLibrary.main(base,combinations)
	if 'cluster' in graphs:
		GraphClusterLibrary.main(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine,methFlank)

def main():
	# Collect arguments
	args = get_args()
	
	# Integer parameters
	num = args.total
	uce = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	methCovThresh = args.thresholdcoverage
	methPerThresh = args.thresholdpercentage
	methFlank = args.methylationflank

	# Arguments for the = bin calibration
	base = args.base
	combinations = args.combinations

	# Element and methylation files
	eFiles = [line.strip() for line in args.efile]
	mFiles = [line.strip() for line in args.mfile]

	# Genome files from UCSC
	sizeGenome = args.genome
	faGenome = args.fasta

	# Lists with the types and directions to use
	typeList = args.elementype
	dirList = args.elementdirection
	nucLine = args.nucleotideline

	# Reverse complement argument
	revCom = args.reversecomplement

	# Which plots to run
	graphs = args.plots

	# A string to add to the out file name in case you want to set up runs and let be
	stringName = args.stringname

	# for each element file provided
	for fileName in eFiles:
		# get the region info to work with
		rangeFeatures = ElementLibrary.main(num,uce,inuce,window,binDir,fileName,sizeGenome,faGenome)
		# separate by direction
		directionFeatures = DirectionLibrary.main(rangeFeatures,fileName,binDir)

		# All elements
		if 'all' in typeList:
			typeList.remove('all')
			allWindow, allNames = FangsLibrary.main(rangeFeatures['combineString'],rangeFeatures['id'],num,uce,inuce,window,nucLine)
			if any(x in graphs for x in ['methylation','cluster']):
				pdMeth = MethylationLibrary.main(mFiles,rangeFeatures,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
			else:
				pdMeth = None
			plotGraphs(pdMeth,allWindow,allNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format('all',uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank)
			if revCom:
				revMeth,revWindow,revNames = RevCompLibrary.main(directionFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				plotGraphs(revMeth,revWindow,revNames,'revComp_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format('all',uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank)

# 		# By Type
		for type in typeList:
			typeBool,typeMeth,typeWindow,typeNames = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
			plotGraphs(typeMeth,typeWindow,typeNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(type,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank)
			if revCom:
				typercMeth,typercWindow,typercNames = RevCompLibrary.main(typeBool,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				plotGraphs(typercMeth,typercWindow,typercNames,'revComp_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(type,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank)

		# By Direction
		for dir in dirList:
			dirBool,dirMeth,dirWindow,dirNames = groupSeparate(dir,directionFeatures,'compareBoundaries',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
			plotGraphs(dirMeth,dirWindow,dirNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}'.format('all',dir,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank)
			for type in typeList:
				typeBool,typeMeth,typeWindow,typeNames = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				plotGraphs(typeMeth,typeWindow,typeNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}'.format(type,dir,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations,methFlank)

if __name__ == "__main__":
	main()