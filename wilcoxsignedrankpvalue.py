"""
Script to get the smallest pval before pval becomes 0

Wren Saylor
March 2018

Copyright 2018 Harvard University, Wu Lab

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
from scipy import stats
import pandas as pd
import numpy as np

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument('-m',"--arraymean",default=1,type=int,help='the largest mean to make')
	parser.add_argument('-n',"--arraynumber",type=int,default=10000,help='the number of elements in the array')
	parser.add_argument('-s',"--arraystep",type=int,default=0.001,help='the size of the steps between elements in the array')
	return parser.parse_args()

# save panda to file
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=False,header=True)

def main():
	args = get_args()
	arraymean = args.arraymean
	arraynumber = args.arraynumber
	arraystep = args.arraystep
	
	collect = []
	
	randomone = stats.norm.rvs(loc=0,scale=1,size=arraynumber)
	for i in np.arange(0,arraymean,arraystep):
		randomtwo = stats.norm.rvs(loc=i,scale=1,size=arraynumber)
		statcoef,statpval = stats.wilcoxon(randomone,randomtwo)
		statpd = pd.DataFrame({'coef':[statcoef],'pval':[statpval],'mean1':[0],'mean2':[i]})
		collect.append(statpd)
	
	concat = pd.concat(collect)
	save_panda(concat,'Wilcox_Signed_Rank_pvalues.txt')

if __name__ == "__main__":
	main()