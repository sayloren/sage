"""
Script to read in tif stack, do some localization counting, return some stats

Wren Saylor
December 2017

Jumana's requests
Short term ideas/ goals: To simply be able to change simple things like:
-As you mentioned file names - CHECK
-number of channels - CHECK
-can compare single copy locus to heterochromatic repeats (intensity: noise threshold will be different, as repeats will always have a higher signal:noise)
-drosophila cells vs. mammalian cells ( as I mentioned Drosophila cells are quite smaller than mammalian cells)

Long-term goals:To properly assay pairing: 1-spot/2-spot assay- Or what I like to call "the ultimate pairing assay"
-volume FISH signal to volume of nucleus
-distance range or 'pairing' distribution between two FISH signals
-distance in Z as well as X-Y
-chance of two different FISH targets to 'pair' within the nucleus.

Jelena's requests
a) during import of images:
-assign channels properly (green is green, red is red, and so on)
-have settings for each of the microscopes that we use that defines how big is the pixel in microns (e.g. on Zeiss, Niko, Vutara)
-input images in batches - CHECK
-have easy to implement standard naming of input files, or no naming at all (like which z-stack and channel image is)

b) output files:
-process 3D distances between spots in a Z-stack (center-to-center)
-output histogram with frequency of cells depending on the micron distance between 2 spots
-have a way to distinguish if the distance between 2 spots was measured in a cell that only has 2 those two spots, or a cell that has 3, 4 or more spots
-measure volumen of spots to see if paired spots are bigger than unpaired

c) stats
-a test to determine if pairing level in one condition is higher than in the other (it could be done with another script, not connected to image processing also)

"""

import argparse
import cv2
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from skimage.morphology import watershed
from scipy import ndimage
# import tiffcapture as tc

# set command line arguments
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("files",type=argparse.FileType('rU'),help="a file with the list of files names to process, only '.' is at the filetype")
	parser.add_argument("-z","--zstackdistance",default="1",type=int,help='distance between the z frames, in micrometers')
	parser.add_argument("-c","--nucleuschannel",default="1",type=str,help='channel number for the nucleus, the other channels will be assumed to be dot counting')
	return parser.parse_args()

# set globals
def set_args(args):
	global filelist
	filelist = [line.strip() for line in args.files]
	global zstackdistance
	zstackdistance = args.zstackdistance
	global dapchannel
	dapchannel = args.nucleuschannel

# splitting the filenames into a dataframe with sample, z stack, and channel in separate columns
def parse_file_names(df):
	pddf = pd.DataFrame(df)
	pddf.columns = ['filename']
	pddf['tempsample'],pddf['filetype'] = pddf['filename'].str.split('.',1).str
	pddf['channel'] = pddf['tempsample'].str.extract('^.*[c|C](.*?$)',expand=True)
	pddf['zstack'] = pddf['tempsample'].str.extract('^.*[z|Z](.*?)[c|C].*?$',expand=True)
	pddf['sample'] = pddf['tempsample'].str.extract('(^.*?)[z|Z].*?$',expand=True)
	outdf = pddf[['filename','sample','zstack','channel']]
	return outdf

# read image
def read_image(filename):
	return cv2.imread(filename)# 0 for grayscale

# convert image to grey scale
def convert_grayscale(image):
	return cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# lighten image with gamma correction
def adjust_gamma(image,gamma):
	#https://www.pyimagesearch.com/2015/10/05/opencv-gamma-correction/
	# build a lookup table mapping the pixel values [0, 255] to their adjusted gamma values
	invGamma = 1.0 / gamma
	table = np.array([((i / 255.0) ** invGamma) * 255
		for i in np.arange(0, 256)]).astype("uint8")
	# apply gamma correction using the lookup table
	return cv2.LUT(image,table)

# blur to denoise
def blur_signal(image,sigma):
	#https://docs.opencv.org/3.0-beta/modules/imgproc/doc/filtering.html
	return cv2.GaussianBlur(image,(5,5),10)

# watershed for segmentation
def watershed_segmentation_of_cells(image,thresh,gray):
	#https://www.pyimagesearch.com/2015/11/02/watershed-opencv/
	# compute the exact Euclidean distance from every binary pixel to the nearest zero pixel, then find peaks in this distance map
	D = ndimage.distance_transform_edt(thresh)
	localMax = peak_local_max(D, indices=False, min_distance=20,labels=thresh)
	# perform a connected component analysis on the local peaks, using 8-connectivity, then appy the Watershed algorithm
	markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0]
	labels = watershed(-D, markers, mask=thresh)
# 	print("[INFO] {} unique segments found".format(len(np.unique(labels)) - 1))
	# loop over the unique labels returned by the Watershed algorithm
	for label in np.unique(labels):
		# if the label is zero, we are examining the 'background' so simply ignore it
		if label == 0:
			continue
		# otherwise, allocate memory for the label region and draw it on the mask
		mask = np.zeros(gray.shape, dtype="uint8")
		mask[labels == label] = 255
		# detect contours in the mask and grab the largest one
		cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
			cv2.CHAIN_APPROX_SIMPLE)[-2]
		c = max(cnts, key=cv2.contourArea)

		# draw a circle enclosing the object
		((x, y), r) = cv2.minEnclosingCircle(c)
		cv2.circle(image, (int(x), int(y)), int(r), (0, 255, 0), 2)
	return image



# horizontally stack before and after images
def stack_before_and_after_process_images(beforeimage,afterimage):
	return np.hstack((beforeimage,afterimage))

# save image to png
def save_image(image,filename):
	cv2.imwrite('{0}.png'.format(filename),image)





# adjust the contract of the image
def histogram_equalization(image):
	return cv2.equalizeHist(image)

# find cell objects
def find_contours(image):
	im2, contours, hierarchy = cv2.findContours(image,cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
	return cv2.drawContours(im2, contours, -1, (0,255,0), 3)

# order of operations to run histogram equalizing
def run_equalize_image(filename):
	image = read_image(filename)
	gray = convert_grayscale(image)
	equalized = histogram_equalization(gray)
	drawing = find_contours(equalized)
	stacked = stack_before_and_after_process_images(equalized,drawing)
	save_image(stacked,'test')

def main():

# Get and set parameters
	args = get_args()
	set_args(args)
	
	# read in files names
	pdfiles = parse_file_names(filelist)
	
	# group by sample
	groupsamples = pdfiles.groupby(['sample'])
	
	# iterate through groups
	for name,group in groupsamples:
	
		# find the dapi channel/ channel where the nucleus is indicated
		dapigroup = group[group['channel'] == dapchannel]
		
		# read in each filename
		for filename in dapigroup['filename']:
# 			run_equalize_image(filename)

			# read in image
			dapiimage = read_image(filename)
			# change to grey instead of with the channel color
			dapigray = convert_grayscale(dapiimage)
			# gamma correction (camera sensitivity to eye sensitivity)
			dapigamma = adjust_gamma(dapigray,10.0)
			# gaussian blur filter to de-noise
			dapiblur = blur_signal(dapigamma,3)
			# equalize histogram
			dapiequalized = histogram_equalization(dapiblur)
			
			
			
			
			dapistacked = stack_before_and_after_process_images(dapigamma,dapiequalized)
			save_image(dapistacked,'test')





	# Equalized hist
# 			img = cv2.imread(filename,0) # 0 for grayscale
# 			equ = cv2.equalizeHist(img)
# 			res = np.hstack((img,equ))
# 			cv2.imwrite('res.png',res)
	
	# Snippets
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

	# Potential Helpful Websites
	# http://mahotas.readthedocs.io/en/latest/
	# https://docs.opencv.org/master/d9/df8/tutorial_root.html
	# https://docs.opencv.org/master/de/d7a/tutorial_table_of_content_core.html
	# https://docs.opencv.org/master/d7/da8/tutorial_table_of_content_imgproc.html
	# https://www.mathworks.com/company/newsletters/articles/the-watershed-transform-strategies-for-image-segmentation.html?refresh=true
	# http://imagej.net/Nuclei_Watershed_Separation
	# https://docs.opencv.org/3.1.0/d5/daf/tutorial_py_histogram_equalization.html
	# https://docs.opencv.org/2.4.13.4/doc/tutorials/imgproc/table_of_content_imgproc/table_of_content_imgproc.html
	# https://mlbernauer.wordpress.com/2014/12/12/statistical-method-for-automated-cell-detection-and-counting/
	# https://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/
	# https://stackoverflow.com/questions/7589012/combining-two-images-with-opencv

if __name__ == "__main__":
	main()
