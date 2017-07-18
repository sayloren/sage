"""
Script to call and execute different sets of analyses for Element Structure

Wren Saylor
July 5 2017

To Do:
Print output graphs to new folder
Have input methylation data in separate callable folder
If no methylation plotting, don't process methylation
Add verbose as a setting to args

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

# set command line arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("mfile", type=argparse.FileType('rU'), help="A file containing a list of paths to the methylation files with unique names separated by newlines, with data for methylation position (chr, start,stop) and methylation % as fourth column'")

	# Genome Files
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
	#parser.add_argument("-n", "--nucleosome", type=str, help="A bedgraph file with data for nucleosome positions, form 'chr, start, stop, occupancy'")
	#parser.add_argument("-s", "--snp", type=str, help="A file with data for snps, form 'chr, start, stop(start+size alt-mutated), ref, ref_size, alt, alt_size, af_adj'")
	parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")

	# Integer Parameters
	parser.add_argument("-t", "--total", type=int, default="600", help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e", "--element", type=int, default="200", help='size of your element (region without flanks), should be an even number')
	parser.add_argument("-i", "--inset", type=int, default="50", help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w", "--window", type=int, default="11", help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b", "--bin", type=int, default="30", help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-mc", "--thresholdcoverage", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')
	parser.add_argument("-mp", "--thresholdpercentage", type=int, default="10", help='size to threshold % methylation data')

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
# 	parser.add_argument('-v',"--verbose",help='print statements')
	return parser.parse_args()

# for type and direction, separate the groups and run the analyses
def groupSeparate(List,directionFeatures,typecolumn,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome):

	# subset by bool presence
	bool = (directionFeatures[directionFeatures[typecolumn] == List])

	# if there is nothing in that set, skip
	if len(bool.index) != 0:
		Meth,Window,Names = TypeLibrary.main(bool,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
	return bool,Meth,Window,Names

# the plotting options, if in the list of plot flags, run graph
def plotGraphs(pdMeth,slidingWinDF,names,fileName,num,uce,inuce,window,graphs,nucLine,base,combinations):
	if 'fangs' in graphs:
		GraphFangLibrary.main(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine)
	if 'signal' in graphs:
		GraphSignalLibrary.main(slidingWinDF,fileName,num,uce,inuce,window,nucLine)
	if 'methylation' in graphs:
		GraphMethLibrary.main(slidingWinDF,pdMeth,fileName,num,uce,inuce,window)
	if 'interactive' in graphs:
		BokehLibrary.main(slidingWinDF,fileName,num,uce,inuce,window,nucLine)
# 	if 'table' in graphs:
# 		GraphTableLibrary.main(slidingWinDF,fileName,num,uce,inuce,window)
# 	if 'combinations' in graphs:
# 		BinLibrary.main(base,combinations)
	if 'cluster' in graphs:
		GraphClusterLibrary.main(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine)

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
			pdMeth = MethylationLibrary.main(mFiles,rangeFeatures,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
			plotGraphs(pdMeth,allWindow,allNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format('all',uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations)
			if revCom:
				revMeth,revWindow,revNames = RevCompLibrary.main(directionFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
				plotGraphs(revMeth,revWindow,revNames,'revComp_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format('all',uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations)

# 		# By Type
		for type in typeList:
			typeBool,typeMeth,typeWindow,typeNames = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
			plotGraphs(typeMeth,typeWindow,typeNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(type,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations)
			if revCom:
				typercMeth,typercWindow,typercNames = RevCompLibrary.main(typeBool,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
				plotGraphs(typercMeth,typercWindow,typercNames,'revComp_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(type,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations)

		# By Direction
		for dir in dirList:
			dirBool,dirMeth,dirWindow,dirNames = groupSeparate(dir,directionFeatures,'compareBoundaries',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
			plotGraphs(dirMeth,dirWindow,dirNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}'.format('all',dir,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations)
			for type in typeList:
				typeBool,typeMeth,typeWindow,typeNames = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome)
				plotGraphs(typeMeth,typeWindow,typeNames,'{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}'.format(type,dir,uce,inuce,num,binDir,window,fileName,stringName),num,uce,inuce,window,graphs,nucLine,base,combinations)

if __name__ == "__main__":
	main()