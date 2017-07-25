"""
Script to graph Signal analyses

Wren Saylor
July 5 2017
"""

import argparse
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from numpy import sin, linspace, pi
import numdifftools as nd
import numpy as np
import pandas as pd
# from peakdetect import peakdetect
import scipy.stats as ss
import scipy as sp
import scipy.fftpack
from scipy.interpolate import splrep, splev
from scipy import signal
from scipy.stats import mstats
import seaborn as sns

# Get just the elemenet
def justElement(region,num,uce,halfwindow,window):
	element = region[(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)]
	return element

# Get just the downstream flank
def downFlank(region,num,uce,halfwindow,window):
	dFlank = region[(((num-uce)/2)+uce-halfwindow):(num-window-window)]
	return dFlank

# Get just the upstream flank
def upFlank(region,num,uce,halfwindow,window):
	uFlank = region[window:(((num-uce)/2)-halfwindow)]
	return uFlank

# Perform fourier transform
def performFourier(region):
	Fs = 1.0 # sampling rate
	Ts = 1.0/Fs # sampling interval
	
	# length of the signal
	nsd = len(region)
	ksd = np.arange(nsd)
	Tsd = nsd/Fs
	
	# two sides frequency range
	frqsd = ksd/Tsd
	
	# one side frequncy range
	frqsd = frqsd[range(nsd/2)]
	
	# fft computing and normalization
	Ysd = np.fft.fft(region)/nsd
	Ysd = Ysd[range(nsd/2)]
	
	return frqsd, Ysd

# Make signal graphs
def graphSignal(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine):
	
	# Parameters used thougout
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)

	# Get mean for AT
	ATNames = [names.index(i) for i in names if 'A' in i or 'T' in i]
	ATDataFrames = [slidingWinDF[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()

	# File name
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)

	# Filename
	pp = PdfPages('Signal_{0}.pdf'.format(fileName))

	gs = gridspec.GridSpec(3,3,height_ratios=[1,1,1])
	gs.update(hspace=.65)

	# Create fitted, first and second derivative lines
	f = splrep(fillX,ATmean,k=5,s=11)
	smoothMean = splev(fillX,f)
	firstDer = splev(fillX,f,der=1)
	firstDer[0:halfwindow] = 0 # small edge effect
	firstDer[-halfwindow:] = 0 # small edge effect
	secondDer = splev(fillX,f,der=2)
	secondDer[0:window] = 0 # small edge effect
	secondDer[-window:] = 0 # small edge effect
	
	#https://blog.ytotech.com/2015/11/01/findpeaks-in-python/
	# return the locations of the inflection points
# 	SDpeaks = peakdetect(secondDer,lookahead=100)
	
	# Plot smoothed mean AT
	ax0 = plt.subplot(gs[0,:])
	ax0.plot(fillX,smoothMean,linewidth=1, color='#3e1638',alpha=0.9)
	ax0.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.set_yticks(ax0.get_yticks()[::2])
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_title('Fitted Mean AT Content',size=8)
	
	# First derivative
	ax5 = plt.subplot(gs[1,:],sharex=ax0)
	ax5.plot(fillX,firstDer,linewidth=1, color='#3e1638',alpha=0.8)
	ax5.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax5.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax5.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax5.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax5.axhline(y=0,linewidth=.1,color='#bd4973',alpha=0.3)
	ax5.set_yticks(ax5.get_yticks()[::2])
	ax5.set_ylabel('Amplitude',size=8)
	ax5.set_title('First Derivative of Fitted Mean',size=8)
	
	# Second derivative
	ax6 = plt.subplot(gs[2,:],sharex=ax0)
	ax6.plot(fillX,secondDer,linewidth=1, color='#3e1638',alpha=0.7)
	#http://stackoverflow.com/questions/13691775/python-pinpointing-the-linear-part-of-a-slope
	#http://stackoverflow.com/questions/16323139/finding-inflection-points-in-spline-fitted-1d-data
	ax6.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax6.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax6.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax6.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax6.axhline(y=0,linewidth=.1,color='#bd4973',alpha=0.3)
	ax6.set_ylabel('Amplitude',size=8)
	ax6.set_xlabel('Position',size=6)
	ax6.set_yticks(ax6.get_yticks()[::2])
	ax6.set_title('Second Derivative of Fitted Mean',size=8)

	sns.despine()
	plt.savefig(pp, format='pdf')
	
	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	
	ax7 = plt.subplot(gs[0,:],sharex=ax0)
	endRange = 100
	widths = np.arange(1, endRange)
	cwtmatr = signal.cwt(firstDer, signal.ricker, widths)
	ax7.imshow(cwtmatr,cmap='RdPu',extent=[0, (num-window), 1, endRange],aspect='auto',vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())
# 	ax7.set_ylabel('Frequency',size=8)
	ax7.set_xlabel('Position',size=6)
	ax7.set_yticks(ax7.get_yticks()[::2])
	ax7.set_title('Continuous Wavelet Transformation Convolved Over Range {0}-{1} for the First Derivative'.format(widths[0],endRange),size=8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	
	gs = gridspec.GridSpec(3,3,height_ratios=[2,1,1])
	gs.update(hspace=.65)
	
	# Short Fourier Transform
	ax8 = plt.subplot(gs[0,:],sharex=ax0)
	sbins = 30
	f1, t1, Zxx1 = signal.stft(firstDer,fs=1.0, window='hann',nperseg=sbins,noverlap=None)#,nperseg=11,noverlap=5
	ax8.pcolormesh(t1,f1,np.abs(Zxx1),cmap='RdPu')
	ax8.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax8.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax8.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#96c85b')
	ax8.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#96c85b')
	ax8.set_ylabel('Frequency',size=8)
	ax8.set_xlabel('Position',size=6)
	ax8.set_yticks(ax8.get_yticks()[::2])
	ax8.set_title('Short Fourier Transform over {0} bins'.format(sbins),size=8)
	
	# First Derivative
	ax9 = plt.subplot(gs[1,:],sharex=ax0)
	ax9.plot(fillX,firstDer,linewidth=1, color='#3e1638')
	ax9.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax9.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax9.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax9.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax9.axvspan((((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow),facecolor = '#ae3e9e',label='',alpha=0.1)
	ax9.axvspan(window,(((num-uce)/2)-halfwindow),facecolor = '#863eae',label='',alpha=0.1)
	ax9.axvspan((((num-uce)/2)+uce-halfwindow),(num-window-window),facecolor = '#ae3e66',label='',alpha=0.1)
	ax9.set_yticks(ax9.get_yticks()[::2])
	ax9.set_xlabel('Position',size=6)
	ax9.set_ylabel('Amplitude',size=8)
	ax9.set_title('First Derivative of Fitted Mean',size=8)
	
	ysdElement = justElement(firstDer,num,uce,halfwindow,window)
	frq2sd, Y2sd = performFourier(ysdElement)
	ysdUp = upFlank(firstDer,num,uce,halfwindow,window)
	frq3sd, Y3sd = performFourier(ysdUp)
	ysdDown = downFlank(firstDer,num,uce,halfwindow,window)
	frq4sd, Y4sd = performFourier(ysdDown)
	
	#FFT for sections of the smoothed second derivative
	ax10 = plt.subplot(gs[2,0])
	ax10.plot(frq3sd,abs(Y3sd),linewidth=1, color='#863eae')
	ax10.set_ylabel('|Y(freq)|',size=8)
	ax10.set_xlabel('Freq(Hz)',size=6) #AT Rate Change
	ax10.set_yticks(ax10.get_yticks()[::2])
	ax11 = plt.subplot(gs[2,1],sharey=ax10)
	plt.setp(ax11.get_yticklabels(), visible=False)
	ax11.plot(frq2sd,abs(Y2sd),linewidth=1, color='#ae3e9e')
	ax11.set_title('Power Series for Highlighted Regions',size=8)# Power Spectrum Analysis for FFT
	ax11.set_xlabel('Freq(Hz)',size=6)
	ax12 = plt.subplot(gs[2,2],sharey=ax10)
	plt.setp(ax12.get_yticklabels(), visible=False)
	ax12.plot(frq4sd,abs(Y4sd),linewidth=1, color='#ae3e66')
	ax12.set_xlabel('Freq(Hz)',size=6)
	
	sns.despine()
	pp.savefig()
	pp.close()

def main(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine):
	graphSignal(slidingWinDF,names,fileName,num,uce,inuce,window,nucLine)
	
if __name__ == "__main__":
	main()