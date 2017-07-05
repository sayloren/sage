"""
Script to call and execute different sets of analyses for Element Structure

Wren Saylor
July 5 2017

To Do:
Print output graphs to new folder
Have input methylation data in separate callable folder
Commenting

Scripts/functions that need cleaning up:
slidingwindow (not have implicit nucleotide combinations)
Graphing (need to have looping for graphs of same style)
Reverse Complement (just messy)
Methylation (some looping could clean things up)

"""

import argparse
import BinLibrary
import BokehLibrary
import DirectionLibrary
import ElementLibrary
import FangsLibrary
import GraphFangLibrary
import GraphMethLibrary
import GraphSignalLibrary
import MethylationLibrary
import RevCompLibrary
import SepDirLibrary
import TypeLibrary

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
	parser.add_argument("-m", "--threshold", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')

	# Directionality Calibration
	parser.add_argument("-n", "--base", type=int, default=2)
	parser.add_argument("-k", "--combinations", type=int, default=30)

	# Specify which groups and graphs to run
	parser.add_argument('-type', "--elementype", default=[], nargs='*', choices=['all','intronic','exonic','intergenic'],help='which group types of element to run') # type
	parser.add_argument('-dir', "--elementdirection", default=[], nargs='*', choices=['+','-','='], help='which group direction of element to run') # direction
	parser.add_argument('-rc', "--reversecomplement",action='store_true', help='if reverse complement sorting required')
	parser.add_argument('-p',"--plots",default=[],nargs='*',choices=['fangs','methylation','signal','interactive'],help='the available graphs to plot')
	return parser.parse_args()

# for type and direction, separate the groups and run the analyses
def groupSeparate(List,directionFeatures,column,fileName,binDir,mFiles,num,uce,inuce,window,methThresh):
	# subset by bool presence
	bool = (directionFeatures[directionFeatures[column] == List])

	# if there is nothing in that set, skip
	if len(bool.index) != 0:
		Meth,Table,Window,CpA,CpT,CpG,CpC,A,T,G,C,Mo = TypeLibrary.main(bool,fileName,binDir,mFiles,num,uce,inuce,window,methThresh)
	return bool,Meth,Table,Window,CpA,CpT,CpG,CpC,A,T,G,C,Mo

# the plotting options, if in the list of plot flags, run graph
def plotGraphs(pdMeth,pdTable,pdWindow,pdCpA,pdCpT,pdCpG,pdCpC,pdA,pdT,pdG,pdC,pdMo,fileName,num,uce,inuce,window,graphs):
	if 'fangs' in graphs:
		GraphFangLibrary.main(pdWindow,pdCpA,pdCpT,pdCpG,pdCpC,pdA,pdT,pdG,pdC,pdMo,fileName,num,uce,inuce,window)
	if 'signal' in graphs:
		GraphSignalLibrary.main(pdWindow,fileName,num,uce,inuce,window)
	if 'methylation' in graphs:
		GraphMethLibrary.main(pdWindow,pdMeth,pdTable,fileName,num,uce,inuce,window)
	if 'interactive' in graphs:
		BokehLibrary.main(pdWindow,fileName,num,uce,inuce,window)

def main():
	# Collect arguments
	args = get_args()
	
	# integer parameters
	num = args.total
	uce = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	methThresh = args.threshold

	# arguments for the = bin calibration
	base = args.base
	combinations = args.combinations

	# element and methylation files
	eFiles = [line.strip() for line in args.efile]
	mFiles = [line.strip() for line in args.mfile]

	# genome files from UCSC
	sizeGenome = args.genome
	faGenome = args.fasta

	# lists with the types and directions to use
	typeList = args.elementype
	dirList = args.elementdirection

	# reverse complement argument
	revCom = args.reversecomplemen

	# which plots to run
	graphs = args.plots

	# for each element file provided
	for fileName in eFiles:
		# get the region info to work with
		rangeFeatures = ElementLibrary.main(num,uce,inuce,window,binDir,fileName,sizeGenome,faGenome)
		# separate by direction
		directionFeatures = DirectionLibrary.main(rangeFeatures,fileName,binDir)

		# All elements
		if 'all' in typeList:
			typeList.remove('all')
			pdWindow,pdCpA,pdCpT,pdCpG,pdCpC,pdA,pdT,pdG,pdC,pdMo = FangsLibrary.main(rangeFeatures,num,uce,inuce,window)
			pdMeth,pdTable = MethylationLibrary.main(mFiles,rangeFeatures,num,uce,inuce,methThresh)
			plotGraphs(pdMeth,pdTable,pdWindow,pdCpA,pdCpT,pdCpG,pdCpC,pdA,pdT,pdG,pdC,pdMo,'{0}_{1}'.format(binDir,fileName),num,uce,inuce,window,graphs)
			# if reverse complement arg is true, run reverse complement
			if revCom:
				pdgroupMeth,pdgroupTable,compWindow,compWinCpA,compWinCpT,compWinCpG,compWinCpC,compWinA,compWinT,compWinG,compWinC,compWinMo = RevCompLibrary.main(directionFeatures,binDir,mFiles,num,uce,inuce,window,methThresh)
				plotGraphs(pdgroupMeth,pdgroupTable,compWindow,compWinCpA,compWinCpT,compWinCpG,compWinCpC,compWinA,compWinT,compWinG,compWinC,compWinMo,'{0}_{1}'.format(binDir,fileName),num,uce,inuce,window,graphs)

# 		BinLibrary.main(base, combinations)

		# By Type
		for type in typeList:
			typeBool,typeMeth,typeTable,typeWindow,typeCpA,typeCpT,typeCpG,typeCpC,typeA,typeT,typeG,typeC,typeMo = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methThresh)
			plotGraphs(typeMeth,typeTable,typeWindow,typeCpA,typeCpT,typeCpG,typeCpC,typeA,typeT,typeG,typeC,typeMo,'{0}_{1}_{2}'.format(type,binDir,fileName),num,uce,inuce,window,graphs)
			if revCom:
				typercMeth,typercTable,typercWindow,typercCpA,typercCpT,typercCpG,typercCpC,typercA,typercT,typercG,typercC,typercMo = RevCompLibrary.main(typeBool,binDir,mFiles,num,uce,inuce,window,methThresh)
				plotGraphs(typercMeth,typercTable,typercWindow,typercCpA,typercCpT,typercCpG,typercCpC,typercA,typercT,typercG,typercC,typercMo,'revComp_{0}_{1}_{2}'.format(type,binDir,fileName),num,uce,inuce,window,graphs)

		# By Direction
		for dir in dirList:
			dirBool,dirMeth,dirTable,dirWindow,dirCpA,dirCpT,dirCpG,dirCpC,dirA,dirT,dirG,dirC,dirMo = groupSeparate(dir,directionFeatures,'compareBoundaries',fileName,binDir,mFiles,num,uce,inuce,window,methThresh)
			plotGraphs(dirMeth,dirTable,dirWindow,dirCpA,dirCpT,dirCpG,dirCpC,dirA,dirT,dirG,dirC,dirMo,'{0}_{1}_{2}'.format(dir,binDir,fileName),num,uce,inuce,window,graphs)

if __name__ == "__main__":
	main()