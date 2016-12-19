"""
Script to make windows using pybedtools, by giving a file with your window, removing those
that contain Ns, so that, nonN regions can be accounted for and
windows containing nonN regions will not have artificially low coverage.
Then calculate coverage of each windows by a feature file.

Wren Saylor, adapeted from Ruth's density script
December 16 2016

Requires bedtools genome file, in same directory as script.

Requires all bed files to have 0-based starts.

"""
import argparse
import pybedtools as pbt
import pandas as pd


def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("-w", "--windows", type=str)
	parser.add_argument("-s", "--size",type=int,help='Size of the intervals in your window file')
	parser.add_argument("file", type=argparse.FileType('rU'),help='A file containing a list'
	'of paths to the feature files you want to process, '
	'separated by newlines')
	parser.add_argument("-n", "--nonN", type=str)
	return parser.parse_args()

def getFileNames(args):
	arFilenames = [line.strip() for line in args.file]
	return arFilenames

def getFeatures(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# original from Ruth's script - the nonNs were used to make windows of certain size
def makeEqualWindows(args):
	btWindows = pbt.BedTool().window_maker(b=str(args.nonN), w = args.windows)
	intWindowLength = args.bin
	btEqualWindows = btWindows.filter(lambda x: len(x) == intWindowLength)
	btSavedEqualWindows = btEqualWindows.saveas('EqualWindowsTemp.bed')
	btEqualWindowsNotStream = getFeatures('EqualWindowsTemp.bed')
	return btEqualWindowsNotStream

def getWindow(args):
	originalWindows = getFeatures(args.windows)
	nonNWindows = getFeatures(args.nonN)
	btWindows = originalWindows.subtract(nonNWindows)
	intWindowLength = args.size
	btEqualWindows = btWindows.filter(lambda x: len(x) == intWindowLength)
	btSavedEqualWindows = btEqualWindows.saveas('EqualWindowsTemp.bed')
	btEqualWindowsNotStream = getFeatures('EqualWindowsTemp.bed')
	return btEqualWindowsNotStream

def coverage(btWindows, btFeatures):
	btCoverage = btWindows.coverage(btFeatures)
	return btCoverage

def cutToFraction(btObject):
	btCutToFraction = btObject.cut([0,1,2,6])
	return btCutToFraction

def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

def getDensity(btWindows, btFeatures):
	btCoverage = coverage(btWindows, btFeatures)
	btCutToFraction = cutToFraction(btCoverage)
	return btCoverage, btCutToFraction

def forEachFeatureFile(strFileName, btWindows):
	btFeatures = getFeatures(strFileName)
	btCoverage, btCutToFraction = getDensity(btWindows, btFeatures)
	saveBedTool(btCoverage, str('Full_coverage_{0}.bed'.format(strFileName)))
	saveBedTool(btCutToFraction, str('Density_for_matrix_{0}.bed'.format(strFileName)))
	return btCoverage, btCutToFraction

def main():
	print 'Getting arguments'
	args = get_args()
	print 'Making windows'
	btWindows = getWindow(args)
	# can use the output bedtool of just the bins to interest with data already have for those bins
	saveBedTool(btWindows, 'Bins_for_density_0based.bed')
	arFilenames = getFileNames(args)
	for filename in arFilenames:
		print 'Getting density files for {0}'.format(filename)
		btCoverage, btCutToFraction = forEachFeatureFile(filename, btWindows)


if __name__ == "__main__":
	main()


