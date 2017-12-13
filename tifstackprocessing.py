"""
Script to read in tif stack, do some localization counting, return some stats

Wren Saylor
December 2017

Jumana's requests
Short term ideas/ goals: To simply be able to change simple things like:
-As you mentioned file names
-number of channels
-can compare single copy locus to heterochromatic repeats (intensity: noise threshold will be different, as repeats will always have a higher signal:noise)
-drosophila cells vs. mammalian cells ( as I mentioned Drosophila cells are quite smaller than mammalian cells)

Long-term goals:To properly assay pairing: 1-spot/2-spot assay- Or what I like to call "the ultimate pairing assay"
-volume FISH signal to volume of nucleus
-distance range or 'pairing' distribution between two FISH signals
-distance in Z as well as X-Y
-chance of two different FISH targets to 'pair' within the nucleus.
"""

import argparse
import cv2
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
# import tiffcapture as tc

# set command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("files",type=argparse.FileType('rU'),help="a file with the list of files names to process, only '.' is at the filetype")
	parser.add_argument("-z", "--zstackdistance",default="1",type=int,help='distance between the z frames, in micrometers')
	return parser.parse_args()

# set globals
def set_args(args):
	global filelist
	filelist = [line.strip() for line in args.files]
	global zstackdistance
	zstackdistance = args.zstackdistance

# splitting the filenames into a dataframe with sample, z stack, and channel in separate columns
def parse_file_names(df):
	pddf = pd.DataFrame(df)
	pddf.columns = ['filename']
	pddf['tempsample'],pddf['filetype'] = pddf['filename'].str.split('.',1).str
	pddf['channel'] = pddf['tempsample'].str.extract('^.*([c|C].*?$)',expand=True)
	pddf['zstack'] = pddf['tempsample'].str.extract('^.*([z|Z].*?)[c|C].*?$',expand=True)
	pddf['sample'] = pddf['tempsample'].str.extract('(^.*?)[z|Z].*?$',expand=True)
	outdf = pddf[['filename','sample','zstack','channel']]
	return outdf

def main():

# Get and set parameters
	args = get_args()
	set_args(args)
	
	# read in files names
	pdfiles = parse_file_names(filelist)
	
	# group by sample
	groupsamples = pdfiles.groupby(['sample'])
	
	# iterate through groups
	# https://docs.opencv.org/3.1.0/d5/daf/tutorial_py_histogram_equalization.html
	# https://docs.opencv.org/2.4.13.4/doc/tutorials/imgproc/table_of_content_imgproc/table_of_content_imgproc.html
	for name, group in groupsamples:
		for filename in group['filename']:
			img = cv2.imread(filename,0) # 0 for grayscale
			equ = cv2.equalizeHist(img)
			res = np.hstack((img,equ))
			cv2.imwrite('res.png',res)
			
			
			
# 			i = cv2.imread(filename)
# 			cv2.imshow(filename,i)
# 			cv2.waitKey(0)
# 			cv2.destroyAllWindows()
			
			# https://pypi.python.org/pypi/TiffCapture
# 			i = tc.opentiff(filename)
# 			_, first_img = i.retrieve()
# 			cv2.namedWindow('video')
# 			for img in i:
# 				temping = cv2.absdiff(first_img, img)
# 				_, temping = cv2.threshold(temping, 5, 255, cv2.THRESH_BINARY)
# 				cv2.imshow('video',temping)
# 				cv2.waitKey(0)
# 			cv2.destroyWindow(name)


	# 2) info 
		# a) inputs args used; pixel conversion and z step size
	# 3) imaging processing (contract, edges, hulls)
		# a) background correction (tiny spots, large spots)
		# b) smoothing
		# c) segmentation
	# 4) locate channel with nucleus, get some info on volume and area
		# a) label and count nuclei, diameter, stats
		# b) measure nuclear shape
	# 5) within nucleus, locate signal in other channels
		# a) foci diameter, foci count, foci distances, centroids, stats
	# 6) print some pictures
	
if __name__ == "__main__":
	main()
