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
import scipy as sp
import scipy.fftpack
from scipy.interpolate import splrep, splev
from scipy import signal as ssignal
import scipy.stats as ss
from scipy.stats import mstats
import seaborn as sns

# Make signal graphs
def graphSignal(slidingWinDF,fileName,num,uce,inuce,window,nucLine):
	
	# Parameters used thougout
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)
	
	# Plot settings
	sns.set_style('ticks')
	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1]) # for the ratios of the graphs against each other
	gs.update(hspace=.8) # setting the space between the graphs
	info = str(fileName) + ', '+ str(len(pdWindow)) + ' - ' "UCES"
	plt.suptitle(info,fontsize=10)

	# Get mean for AT
	ATNames = [names.index(i) for i in names if 'A' in i or 'T' in i]
	ATDataFrames = [slidingWinDF[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()

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
	
	# Plot smoothed mean AT
	ax4 = plt.subplot(gs[0,:])
	ax4.plot(fillX,smoothMean,linewidth=1, color='#3e1638',alpha=0.9)
	ax4.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax4.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax4.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax4.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax4.set_yticks(ax4.get_yticks()[::2])
	ax4.set_ylabel('% AT Content',size=8)
	ax4.set_title('Fitted Mean AT Content',size=8)
	
	# First derivative
	ax5 = plt.subplot(gs[1,:],sharex=ax4)
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
	ax6 = plt.subplot(gs[2,:],sharex=ax4)
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
	
	gs = gridspec.GridSpec(3,3,height_ratios=[2,1,1])
	gs.update(hspace=.65)
	
	# Short Fourier Transform
	ax7 = plt.subplot(gs[0,:],sharex=ax4)
	f1, t1, Zxx1 = ssignal.stft(firstDer,fs=1.0, window='hann',nperseg=30,noverlap=None)#,nperseg=11,noverlap=5
	ax7.pcolormesh(t1,f1,np.abs(Zxx1),cmap='RdPu')
	ax7.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax7.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax7.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#96c85b')
	ax7.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#96c85b')
	ax7.set_ylabel('Frequency',size=8)
	ax7.set_xlabel('Position',size=6)
	ax7.set_yticks(ax7.get_yticks()[::2])
	ax7.set_title('Short Fourier Transform',size=8)#30 bp bins
	
	# First Derivative
	ax8 = plt.subplot(gs[1,:],sharex=ax4)
	ax8.plot(fillX,firstDer,linewidth=1, color='#3e1638')
	ax8.axvline(x=(((num-uce)/2)+(inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax8.axvline(x=(((num-uce)/2)+(uce-inuce-halfwindow)),linewidth=.05,linestyle='dashed',color='#e7298a')
	ax8.axvline(x=(((num-uce)/2)-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax8.axvline(x=(((num-uce)/2)+uce-halfwindow),linewidth=.05,linestyle='dashed',color='#bd4973')
	ax8.axvspan((((num-uce)/2)-halfwindow),(((num-uce)/2)+uce-halfwindow),facecolor = '#ae3e9e',label='',alpha=0.1)
	ax8.axvspan(window,(((num-uce)/2)-halfwindow),facecolor = '#863eae',label='',alpha=0.1)
	ax8.axvspan((((num-uce)/2)+uce-halfwindow),(num-window-window),facecolor = '#ae3e66',label='',alpha=0.1)
	ax8.set_yticks(ax8.get_yticks()[::2])
	ax8.set_xlabel('Position',size=6)
	ax8.set_ylabel('Amplitude',size=8)
	ax8.set_title('First Derivative of Fitted Mean',size=8)
	
	Fs = 1.0 # sampling rate
	Ts = 1.0/Fs # sampling interval
	y2sd = firstDer[(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)]
	n2sd = len(y2sd) # length of the signal
	k2sd = np.arange(n2sd)
	T2sd = n2sd/Fs
	frq2sd = k2sd/T2sd # two sides frequency range
	frq2sd = frq2sd[range(n2sd/2)] # one side frequncy range
	Y2sd = np.fft.fft(y2sd)/n2sd # fft computing and normalization
	Y2sd = Y2sd[range(n2sd/2)]
	y3sd = firstDer[window:(((num-uce)/2)-halfwindow)]
	n3sd = len(y3sd)
	k3sd = np.arange(n3sd)
	T3sd = n3sd/Fs
	frq3sd = k3sd/T3sd
	frq3sd = frq3sd[range(n3sd/2)]
	Y3sd = np.fft.fft(y3sd)/n3sd
	Y3sd = Y3sd[range(n3sd/2)]
	y4sd = firstDer[(((num-uce)/2)+uce-halfwindow):(num-window-window)]
	n4sd = len(y4sd)
	k4sd = np.arange(n4sd)
	T4sd = n4sd/Fs
	frq4sd = k4sd/T4sd
	frq4sd = frq4sd[range(n4sd/2)]
	Y4sd = np.fft.fft(y4sd)/n4sd
	Y4sd = Y4sd[range(n4sd/2)]
	
# 	FFT for sections of the smoothed second derivative
	ax9 = plt.subplot(gs[2,0])
	ax9.plot(frq3sd,abs(Y3sd),linewidth=1, color='#863eae')
	ax9.set_ylabel('|Y(freq)|',size=8)
	ax9.set_xlabel('Freq(Hz)',size=6)#AT Rate Change
	ax9.set_yticks(ax9.get_yticks()[::2])
	ax10 = plt.subplot(gs[2,1],sharey=ax9)
	plt.setp(ax10.get_yticklabels(), visible=False)
	ax10.plot(frq2sd,abs(Y2sd),linewidth=1, color='#ae3e9e')
	ax10.set_title('Power Series for Highlighted Regions',size=8)# Power Spectrum Analysis for FFT
	ax10.set_xlabel('Freq(Hz)',size=6)
	ax11 = plt.subplot(gs[2,2],sharey=ax9)
	plt.setp(ax11.get_yticklabels(), visible=False)
	ax11.plot(frq4sd,abs(Y4sd),linewidth=1, color='#ae3e66')
	ax11.set_xlabel('Freq(Hz)',size=6)
	sns.despine()
	pp.savefig()
	
	sns.despine()
	pp.savefig()
	pp.close()

def main(slidingWinDF,fileName,num,uce,inuce,window,nucLine):
	graphSignal(slidingWinDF,fileName,num,uce,inuce,window,nucLine)
	
if __name__ == "__main__":
	main()