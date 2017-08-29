"""
Script to take in a list of element files, and separate their regions into genic, exonic, intronic, and intergenic
Uses a good deal of Ruth's Create_boundaries script

Wren Saylor
August 29 2017
"""

import argparse
import pybedtools as pbt

def get_args():
	# Catch command-line arguments
	parser = argparse.ArgumentParser(description="Script to return interval files for exons, introns, intergenic, genic regions")
	parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files')
	parser.add_argument('-n', '--nonN', type=str,help="interval file for nonN regions, 0-based starts")
	parser.add_argument('-e', '--exons', type=str,help="exonic intervals filename, 0-based starts")
	parser.add_argument('-i', '--introns', type=str,help="introns intervals filename, 0-based stars")
	return parser.parse_args()

def getFeatures(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

def saveBedTool(btObject, strFilename):
	btObject.saveas(strFilename)

def sortAndMerge(btObject):
	btSorted = btObject.sort()
	btMerged = btSorted.merge()
	return btMerged

def main():
	args = get_args()

	#Obtain bedtool objects
	btExons = getFeatures(str(args.exons))
	btIntrons = getFeatures(str(args.introns))
	btnonN = getFeatures(str(args.nonN))

	#Merge (also called collapse) the exon, intron, and nonN bedfiles
	print "Merging exons"
	btExonsFlat = sortAndMerge(btExons)

	print "Merging introns"
	btIntronsFlat = sortAndMerge(btIntrons)

	btnonNFlat = sortAndMerge(btnonN)

	#Intersect the exons and introns with the nonN regions to produce only regions with no N's (you would expect there
	#would be no exons or introns in nonN regions, but just to remove any chromosomes we're not interested in like chrM
	btExonsNonN = btExonsFlat.intersect(btnonNFlat)
	btIntronsNonN = btIntronsFlat.intersect(btnonNFlat)

	#Subtract the exons from the introns
	btIntronsOnly = btIntronsNonN.subtract(btExonsNonN)
	btIntronsOnlyFlat = sortAndMerge(btIntronsOnly)

	#Concatenate and then merge exons and introns
	btGenicnonN =btExonsNonN.cat(btIntronsOnlyFlat)
	btGenicFlat = btGenicnonN.merge()

	#Subtract genic regions from full nonN file to give intergenic regions
	btIntergenicnonN = btnonN.subtract(btGenicFlat)
	btIntergenicFlat = sortAndMerge(btIntergenicnonN)

	eFiles = [line.strip() for line in args.file]
	for fileName in eFiles:
		btFeatures = getFeatures(fileName)
		genicFeatures = btFeatures.intersect(btGenicFlat)
		exonicFeatures = btFeatures.intersect(btExonsNonN)
		intronicFeatures = genicFeatures.subtract(exonicFeatures)
		intergenicFeatures = genicFeatures.subtract(btFeatures)
		#Output these 0-based files
		saveBedTool(exonicFeatures, 'Exons_0based_{0}.bed'.format(fileName))
		saveBedTool(intronicFeatures, 'Introns_0based_{0}.bed'.format(fileName))
		saveBedTool(genicFeatures, 'Genic_0based_{0}.bed'.format(fileName))
		saveBedTool(intergenicFeatures, 'Intergenic_0based_{0}.bed'.format(fileName))


if __name__ == "__main__":
	main()

