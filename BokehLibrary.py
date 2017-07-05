"""

Script to run the same analysises as Fang element library, but makes bokeh plots

Wren Saylor
June 20 2017

"""

import argparse
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, curdoc
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Select
from bokeh.models.widgets import Toggle
from numpy import sin, linspace, pi
import pandas as pd
import seaborn as sns

def bokehOut(pdWindow,fileName,num,uce,inuce,window):
	fillX = range(0,(num-window))
	source = ColumnDataSource(data=dict(x=fillX, mean=pdWindow.mean(axis=0), std=pdWindow.std(axis=0)))
	output_file('Fangs_{0}.html'.format(fileName))
	p = figure(plot_width=1500, plot_height=600, min_border=10, min_border_left=50,toolbar_location="above",title="Mean AT Content Across Base Pair Position")
	p.line('x','mean',line_width=2,color='#3e1638',source=source)
	p.yaxis.axis_label = "% AT Content"
	p.xaxis.axis_label = "Nucleotide Postion"
	p.background_fill_color = "#fafafa"
	
	select = Select(title="Option:", value="All", option=["All", "Exonic", "Intronic", "Intergenic"])
	#toggle = Toggle(label="Reverse Complement")
	
	sd = figure(plot_width=1500, plot_height=200, x_range=p.x_range, min_border=10, min_border_left=50,title="Standard Deviation")
	sd.line('x','std',line_width=2,color='#3e1638',alpha=.5,source=source)
	sd.background_fill_color = "#fafafa"
	
	widgets = row(select,toggle)
	show(column(widgets,p,sd))

def main(pdWindow,fileName,num,uce,inuce,window):
	bokehOut(pdWindow,fileName,num,uce,inuce,window)

if __name__ == "__main__":
	main()
