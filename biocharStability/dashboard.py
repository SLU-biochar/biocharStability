"""
-*- coding: utf-8 -*-

biochar stability / dashboard.py

set of utility functions to create interactive bokeh dashboard from the data
"""
import numpy as np
import pandas as pd
from bokeh.plotting import figure, show, save, output_file
from bokeh.resources import CDN
from bokeh.embed import file_html, json_item, components
from bokeh.io import output_notebook 
from bokeh.layouts import column, grid, layout, gridplot, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, DataTable, DateFormatter,
                          TableColumn, LegendItem, Legend,  MultiChoice, Button, DataRange1d, Plot,
                          BasicTicker, Circle, LinearAxis, PanTool,WheelZoomTool, HoverTool, Grid)
from bokeh.palettes import all_palettes

from .utils import select_mineralization_data

def generate_bokeh_app(
    component_ids=[],
    component_js_names=[],
    bokeh_elements=[],
    md_template_file = 'dashboard-app-template.md',
    md_outpout_file = 'bokeh-app/app.md',
    js_output_file_dir = 'bokeh-app/'
):
    '''
    BOKEH - Generate and saves the bokeh components (as javascript files) and the static app (as markdown for hugo website, from template file), based on the list of components passed as argument.

    USAGE: \n
        generate_bokeh_app(
            component_ids=[1, 2], # must be present in the template .md file used as input
            component_js_names=['component_1', 'component_2'], # free to choose the name of the js file
            bokeh_elements=[component_1, component_2], # bokeh objects available in the notebook
            md_outpout_file = 'C:/github/biochar-systems-dev/content/en/stability/app.md',
            js_output_file_dir = 'C:/github/biochar-systems-dev/content/en/stability/',
        )

    '''
    with open(md_template_file, 'r') as mdf:
        contents = mdf.readlines()
        # type(contents), len(contents) # list with 1 element for each line of the document
    for component_id, component_js_name, bokeh_elem in zip(component_ids, component_js_names, bokeh_elements):
        script, div = components(bokeh_elem, wrap_script=False)
        with open(js_output_file_dir+component_js_name+'.js', 'w') as f:
            f.write(script) # .split('\n', 1)[1].rsplit('\n', 2)[0] # equivalent to wrap_script=False

        # find the marker in the list
        i=0
        while i < len(contents):
            if 'component_'+str(component_id) in contents[i]:
                print(i, contents[i])
                index=i
                break
            i+=1
        #index = contents.index('<div id="component_'+str(component_id)+'" style="height:90vh;" >\n')
        value = "    "+div
        contents.insert(index+1, value)

        # find the marker in the list
        index = contents.index('<div id="bokeh_component_scripts">\n')
        value = '<script type="text/javascript" src="../'+component_js_name+'.js"></script>\n'
        contents.insert(index+1, value)
        print(">>> bokeh component was inserted here")

    with open(md_outpout_file, "w") as f:
        contents = "".join(contents)
        f.write(contents)


def merge_data_metadata(data, metadata, asCDS=False):
    '''
    Merges the data and metadata as a single dataframe, with lists as cell concent for the timeseries in data.
    
    Returns either a dataframe (if asCDS = False), or a ColumnDataSource (if asCDR=True) useful for bokeh plots.

    USAGE: \n
        datametadata = merge_data_metadata(data, metadata, asCDS=False)
    
    OR \n
        TheSource = merge_data_metadata(data, metadata, asCDS=True)
    '''
    datametadata = metadata.copy(deep=True)
    newcols = data.columns.drop(['ID_obs', 'ID_art', 'date'])
    for c in newcols:
        datametadata[c]=np.zeros(len(datametadata))
        L = []
        for i,_ in datametadata.iterrows():
            L.append(select_mineralization_data(i, data)[c].values)
        
        datametadata[c] = L
    TheSource = ColumnDataSource(datametadata)
    return TheSource if asCDS else datametadata


def make_scatter(xname, yname, xRanges, yRanges, src, x_axis_label='X Axis', y_axis_label='Y Axis', xax=False, yax=False, ):
    '''
    BOKEK - create and return a scatter plot object, which can be inserted in a grid layout. Used for making pairplots (with scatter and histograms).
    
    '''
    mbl = 40 if yax else 0
    mbb = 40 if xax else 0
    plot = figure(
        x_range=xRanges[xname], y_range=yRanges[yname], background_fill_color="#efe8e2",
        x_axis_label=x_axis_label, y_axis_label=y_axis_label, toolbar_location=None,
        border_fill_color='white', width=500 + mbl, height=500 + mbb,
        min_border_left=2+mbl, min_border_right=2, min_border_top=2, min_border_bottom=2+mbb)

    circle = Circle(x=xname, y=yname, fill_color="color", fill_alpha=0.2, size=4, line_color="color")
    r = plot.add_glyph(src, circle)

    xRanges[xname].renderers.append(r)
    yRanges[yname].renderers.append(r)

    xticker = BasicTicker()
    if xax:
        xaxis = LinearAxis()
        xaxis.axis_label = xname
        plot.add_layout(xaxis, 'below')
        xticker = xaxis.ticker
    plot.add_layout(Grid(dimension=0, ticker=xticker))

    yticker = BasicTicker()
    if yax:
        yaxis = LinearAxis()
        yaxis.axis_label = yname
        yaxis.major_label_orientation = 'vertical'
        plot.add_layout(yaxis, 'left')
        yticker = yaxis.ticker
    plot.add_layout(Grid(dimension=1, ticker=yticker))

    plot.add_tools(PanTool(), WheelZoomTool())
    
    hover = HoverTool(tooltips = [('Observation', '@ID_obs'),
                                 ])
    plot.add_tools(hover)

    return plot


def make_histogram(x, df, x_axis_label='X Axis', y_axis_label='Y Axis', xax=False, yax=False):
    '''
    BOKEH - create and return a scatter histogram object, where histogram bins are built with numpy
    '''
    # calculate histogram properties using numpy
    arr_hist, edges = np.histogram(df[x].dropna(), bins = int(np.sqrt(len(df[x].dropna()))) )
    # Put the information in a dataframe
    dd = pd.DataFrame({'tops': arr_hist, 
                       'left': edges[:-1], 
                       'right': edges[1:]})
    dd['interval'] = ['%d to %d' % (left, right) for left, right in zip(dd['left'], dd['right'])]
    src = ColumnDataSource(dd)
    
    # Create the blank plot
    p = figure(plot_height = 600, plot_width = 600, 
              x_axis_label = x_axis_label, 
               y_axis_label = y_axis_label)
    # Add a quad glyph with source this time
    p.quad(bottom=0, top='tops', left='left', right='right', source=src,
           fill_color='red', line_color='black', fill_alpha = 0.75,
           hover_fill_alpha = 1.0, hover_fill_color = 'navy')
    # Add a hover tool referring to the formatted columns
    hover = HoverTool(tooltips = [(x_axis_label, '@interval'),
                                 ('Count', '@tops')])
    # Style the plot
    #p = style(p)
    # Add the hover tool to the graph
    p.add_tools(hover)

    # Show the plot
    return p

