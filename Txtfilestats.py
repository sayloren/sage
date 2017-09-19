"""
Script to process the sentence structure in a .txt file (or some other simple plain text file type)

Wren Saylor
September 18 2017

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
import pandas as pd
import re
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=argparse.FileType('rU'),help='A file containing a list of paths to the files with unique names separated by newlines')
	# total, or by paragraph https://stackoverflow.com/questions/44288169/parse-text-file-python-and-covert-to-pandas-dataframe
	parser.add_argument('-s',"--searchstring",type=str,help='if there is a certain string you would like to count, besides the number of commas')
	parser.add_argument("-n","--numberlongest",type=int,default="3",help='the number of longests sentences to return')
	parser.add_argument('-i', "--ignorebrackets",action='store_true',help='it you want to ignore information between brackets')
	return parser.parse_args()

# open and read in file
def open_and_read(file):
	openfile = open(file)
	readfile = openfile.read()
	return readfile

# regex {}
def remove_inside_curly_bracket(text):
	nocurlybracket = re.sub('{.*?}', '', text)
	return nocurlybracket

# break by sentence
def split_into_sentences(text):
# from https://stackoverflow.com/questions/4576077/python-split-text-on-sentences

	# declaration of regular expressions
	caps = "([A-Z])"
	prefixes = "(Mr|St|Mrs|Ms|Dr)[.]"
	suffixes = "(Inc|Ltd|Jr|Sr|Co)"
	starters = "(Mr|Mrs|Ms|Dr|He\s|She\s|It\s|They\s|Their\s|Our\s|We\s|But\s|However\s|That\s|This\s|Wherever)"
	acronyms = "([A-Z][.][A-Z][.](?:[A-Z][.])?)"
	websites = "[.](com|net|org|io|gov)"
	digits = "([0-9])"

	text = " " + text + "  "
	text = text.replace("\n"," ")
	# does not split the line at decimals such as 5.5
	text = re.sub(digits + "[.]" + digits,"\\1<prd>\\2",text)
	# if "e.g." in text
	text = text.replace("e.g.","e<prd>g<prd>")
	# if "i.e." in text
	text = text.replace("i.e.","i<prd>e<prd>")
	# if "et al." in text
	text = text.replace("al.","al<prd>")
	# if a webaddress is present between ()
	text = re.sub(r'\(https?:\/\/(.*?\))','<webadress>',text,flags=re.MULTILINE) # webaddress within ()
	text = re.sub(prefixes,"\\1<prd>",text)
	text = re.sub(websites,"<prd>\\1",text)
	if "Ph.D" in text: text = text.replace("Ph.D.","Ph<prd>D<prd>")
	text = re.sub("\s" + caps + "[.] "," \\1<prd> ",text)
	text = re.sub(acronyms+" "+starters,"\\1<stop> \\2",text)
	text = re.sub(caps + "[.]" + caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>\\3<prd>",text)
	text = re.sub(caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>",text)
	text = re.sub(" "+suffixes+"[.] "+starters," \\1<stop> \\2",text)
	text = re.sub(" "+suffixes+"[.]"," \\1<prd>",text)
	text = re.sub(" " + caps + "[.]"," \\1<prd>",text)
# 	if """ in text: text = text.replace("."","".")
	if "\"" in text: text = text.replace(".\"","\".")
	if "!" in text: text = text.replace("!\"","\"!")
	if "?" in text: text = text.replace("?\"","\"?")
	text = text.replace(".",".<stop>")
	text = text.replace("?","?<stop>")
	text = text.replace("!","!<stop>")
	text = text.replace("<prd>",".")
	sentences = text.split("<stop>")
	sentences = sentences[:-1]
	sentences = [s.strip() for s in sentences]
	return sentences

# to panda
def list_of_sentences_to_panda(sentencelist,searchstring):
	pddataframe = pd.DataFrame(sentencelist)
	pddataframe.columns = ['sentences']
	pddataframe['comma'] = pddataframe['sentences'].apply(lambda x: str.count(x,','))
	pddataframe['length'] = pddataframe['sentences'].str.split().apply(len)
	if searchstring:
		pddataframe[searchstring] = pddataframe['sentences'].apply(lambda x: str.count(x,searchstring))
	return pddataframe

# get stats
def get_mean_and_sd_from_df(pddataframe,column):
	pdmean = pddataframe[column].mean()
	pdsd = pddataframe[column].std()
	return pdmean, pdsd

# return n longest sentences
def get_longest_sentences(sentencedataframe,numberlongest):
	largestsentences = sentencedataframe.nlargest(numberlongest,'length')
	return largestsentences

# save to panda
def save_panda_to_file(statdf,stringfilename,useindex):
	statdf.to_csv(stringfilename,sep="\t",index=useindex)

# histogram
def graph_hist(pdsentences,column,file):
	sns.set_style('ticks')
	pp = PdfPages('histogram_for_{0}_column_{1}.pdf'.format(file,column))
	sns.set_palette("husl")
	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	gs.update(hspace=.5)
	ax0 = plt.subplot(gs[0])
	sns.distplot(pdsentences[column],ax=ax0)
	plt.title('Sentence Length for {0}'.format(file),size=10)
	sns.set_context(font_scale=.5)
	sns.despine()
# 	plt.tight_layout()
	plt.savefig(pp, format='pdf',bbox_inches='tight')#
	pp.close()

def main():
	# get args
	args = get_args()
	numberlongest = args.numberlongest
	ignorebrackets = args.ignorebrackets
	searchstring = args.searchstring
	infile = [line.strip() for line in args.file]
	
	# for each file
	for file in infile:
		readfile = open_and_read(file)
		
		# if want to ignore data inside curly brackets
		if ignorebrackets:
			readfile = remove_inside_curly_bracket(readfile)

		# split up the sentences
		splitsentences = split_into_sentences(readfile)
		
		# into pd dataframe
		pdsentences = list_of_sentences_to_panda(splitsentences,searchstring)
		
		# get the mean and sd
		lengthmean, lengthsd = get_mean_and_sd_from_df(pdsentences,'length')
		commamean, commasd = get_mean_and_sd_from_df(pdsentences,'comma')
		
		# longest sentences
		largestsentences = get_longest_sentences(pdsentences,numberlongest)
		
		# get stats to df
		sentencestatistics = pd.DataFrame({'mean':[lengthmean,commamean],'sd':[lengthsd,commasd]},index=['length','commas'])
		if searchstring:
			strmean, strsd = get_mean_and_sd_from_df(pdsentences,searchstring)
			searchdf = pd.DataFrame({'mean':[strmean],'sd':[strsd]},index=[searchstring])
			sentencestatistics = sentencestatistics.append(searchdf)
		
		# to output files
		save_panda_to_file(sentencestatistics,'summary_statistics_for_{0}'.format(file),True)
		save_panda_to_file(largestsentences['sentences'],'longest_sentences_for_{0}'.format(file),False)
		save_panda_to_file(pdsentences['sentences'],'check_that_sentences_make_sense_for_{0}'.format(file),True)
		
		# graph
		graph_hist(pdsentences,'length',file)

if __name__ == "__main__":
	main()