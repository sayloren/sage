"""
    
    Script for converting density data from UCSC (was wig, used wig2bed, now bedfile where
    fourth column is density) to the bin-ed density needed to run in Ruth's Density script
    
    Input: A master text file with the file names of the interval files you wish to convert
    and the bin file returned from Ruth's Density script with the windows to use
    
    Output: Density files with proper windows binned that can be piped right into Ruth's Matrix
    maker
    
    To Do: The merge might be over merging the lines - compressing instead of bin-ing
    
    Version 2: added sort
    
    Version 3: streamlined function
    
    Wren Saylor, adapted from Ruth's Density Script
    December 13 2016
    
    """

import argparse
import pybedtools as pbt
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Description")
        parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a'
                            'list of paths to the feature files you want to process, '
                            'separated by newlines')
                            parser.add_argument("-b", "--bin", type=str, default="Bins_for_density_0based.bed"
                                                parser.add_argument("-m", "--mode", type=str, default="mean"))
                            
        return parser.parse_args()

# 1 - read in files
def getFileNames(args):
    arFilenames = [line.strip() for line in args.file]
        return arFilenames

# 2 - get features of file
def getFeatures(strFileName):
    preSortbtFeatures = pbt.BedTool(strFileName)
        btFeatures = preSortbtFeatures.sort()
        return btFeatures

# 3 - left outer join intersect
def intersectWindow(btWindows, btFeatures):
    btIntersect = btWindows.intersect(btFeatures,loj=True)
        return btIntersect

# 4 - grab just columns for window coordinates and densities
def cutToFraction(btIntersect):
    btCutToFraction = btIntersect.cut([0,1,2,6])
        return btCutToFraction

# 5 - merge the density features in the windows by mean
def mergeWindow(btCutToFraction,mode):
    btMerge = btCutToFraction.merge(c=4,o=mode)
        return btMerge

# 6 - save to bedtool
def saveBedTool(btObject, strFilename):
    btObject.saveas(strFilename)

def main():
    args = get_args()
        btWindows = getFeatures(args.bin)
        arFilenames = getFileNames(args)
        mode = args.mode
        for filename in arFilenames:
            btFeatures = getFeatures(filename)
                btIntersect = intersectWindow(btWindows,btFeatures)
                btCutToFraction = cutToFraction(btIntersect)
                btMerge = mergeWindow(btCutToFraction,mode)
                saveBedTool(btMerge, str('Density_for_matrix_{0}.bed'.format(filename)))

if __name__ == "__main__":
    main()

