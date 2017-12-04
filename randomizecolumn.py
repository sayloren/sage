"""
Script to randomize the fourth column of a bed file

Wren Saylor
Deceber 4 2017

Copyright 2017 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

import argparse
import pybedtools as pbt
import pandas as pd
import numpy as np
import random

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing a'
		'list of paths to the feature files you want to process,'
		'separated by newlines')
	parser.add_argument("-c","--chunk",action='store_true',help='if chunk randomize is required, else shift will be used')
	parser.add_argument("-n","--numberoperation",type=int,default="1",help='size of the random movement'
		'to make,suggested options [1,3000,20000,100000,200000]')
	parser.add_argument("-l","--column",type=int,help='column to randomize, 0 based')
	return parser.parse_args()

# 1 - read in files
def get_file_names(args):
	arfilenames = [line.strip() for line in args.file]
	return arfilenames

# 2 - get features of file
def get_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# 3 - bedtool to panda for merge to be implemented correctly - depreciated; used in conjunction with intersect
def bedtoolToPanda(btobject,column):
	pdObject = pd.read_table(btobject.fn, header=None)
	pdObject['density'] = pdObject.loc[:,column]
	pdObject['chr'] = pdObject.loc[:,0]
	pdObject['start'] = pdObject.loc[:,1]
	pdObject['stop'] = pdObject.loc[:,2]
	return pdObject

# 4a - chunk and randomize fourth column if not shift
def chunk_randomize_fourth(btfeatures,introll):
	length = len(btfeatures.index)
	array = btfeatures['density'].values
	chunks = [array[i:i+introll] for i  in range(0, length, introll)] # but this isn't randomly split - problem?
	shuffle = random.sample(chunks, k=len(chunks))
	collect = np.concatenate(shuffle,axis=0)
	btfeatures['chunk'] = collect
	dropfeatures = btfeatures[['chr','start','stop','chunk']]
	return dropfeatures

# 4b - roll the fourth column if shift randomize
def numpy_roll_fourth(btfeatures,introll):
	btfeatures['shift'] = np.roll(btfeatures.density,introll)
	dropfeatures = btfeatures[['chr','start','stop','shift']]
	return dropfeatures

# 5 - save file from panda
def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', header=False, index=False)

def main():
	args = get_args()
	arfilenames = get_file_names(args)
	for filename in arfilenames:
		btfeatures = get_features(filename)
		pandafeatures = bedtoolToPanda(btfeatures,args.column)
		if args.chunk == True:
			randomfeatures = chunk_randomize_fourth(pandafeatures,args.numberoperation)
			randommethod = 'chunk'
		else: 
			randomfeatures = numpy_roll_fourth(pandafeatures,args.numberoperation)
			randommethod = 'shift'
		savePanda(randomfeatures,'{0}_randomize_{1}_{2}.bed'.format(filename,randommethod,args.numberoperation))

if __name__ == "__main__":
	main()