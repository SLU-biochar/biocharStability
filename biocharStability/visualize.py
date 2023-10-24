"""
-*- coding: utf-8 -*-

biochar stability / visualize.py 

set of utility functions to visualise the data and export static figures (matplotlib, ternary)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import seaborn as sns
import ternary as ter
from adjustText import adjust_text
from inspect import getfullargspec, signature
from ast import literal_eval
import re

from .utils import select_mineralization_data, select_timeseries
from .analyse import singleExp, singleExp_u, doubleExp, doubleExp_u, tripleExp, tripleExp_u, powerModel, powerModel_q10, std_fmodel

map_model_function = {
    "singleExp":singleExp,
    "singleExp_u":singleExp_u,
    "doubleExp":doubleExp,
    "doubleExp_u":doubleExp_u,
    "tripleExp":tripleExp,
    "tripleExp_u":tripleExp_u,
    "powerModel":powerModel,
    "k_singleExp":singleExp,
    "k_singleExp_u":singleExp_u,
    "k_doubleExp":doubleExp,
    "k_doubleExp_u":doubleExp_u,
    "k_tripleExp":tripleExp,
    "k_tripleExp_u":tripleExp_u,
    "k_powerModel":powerModel,
}

sns.set_theme()

def make_corr_scatter(df, x='HHT', y='BC_100', c='BiomassClass', m = 'PyrolysisClass', 
                      an=False, fi=None, 
                      fig=None, ax=None, saveFig=False, figsize=(6,6)):
    '''
    Plot of scatter between two columns of the passed dataframe, with options to set color and marker via other columns, and annotate each datapoint.
    ''' 
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    sns.scatterplot(data=df, x=x, y=y,
                    style=m, hue=c, palette=sns.color_palette("Set2", as_cmap=False)[:len(df[c].unique())], 
                    ax = ax, s=80)
    ax.set_xlabel(x, fontsize=14)
    ax.set_ylabel(y, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(facecolor='white', framealpha=0.2, fontsize=8) 
    
    if an is not False:
        texts = an
        xss = list(df[x])
        yss = list(df[y])
        for i, txt in enumerate(texts):
            ax.annotate(txt, xy=(xss[i], yss[i]), xycoords="data", fontsize=6)
    
    if saveFig:
        x = x.replace("/", "")
        fig.savefig('simulations/scatter-'+y+'---'+x+'.png', dpi=300, bbox_inches='tight')

def plot_timeseries_by_group(
    data,
    metadata,
    obs_sets,
    col_to_plot='C_bc_rem_rel',
    exclude_obs = [49, 56],
    titles=['a. Incubations longer than 30 months (17 obs)',
           'b. Incubations between 20 and 30 months (38 obs)',
           'c. Incubations between 12 and 20 months (15 obs)',
           'd. Incubations shorter than 12 months (58 obs)',
         ],
    ylabel='Timeseries', legend_title='Categories',
    factor=1,
    marker_list = ['x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x',]*10,
    color_list = sns.color_palette("hls", 8),
    line=False, lw=0.35,
    nrows=2, ncols=2, figsize=(16,10),
    saveFig=False, pathFig= 'img/timeseries_for_.png', dpiFig=300,
):
    '''
    Using matplotlib, create scatter plots of timeseries from the incubation data, grouped according to the sets of observations passed as argument (e.g. by experimental duration, by type of biomass).

    Note: if timeseries values shall be converted (e.g. decay rates from fraction/day to %/year) perform the conversion on the data. Alternatively use the parameter `factor`.

    Inputs:
    - data: data df
    - metadata: metadata df
    - obs_sets: list of sets of integers, corresponding ot ID_obs 
    - col_to_plot: str, corresponding to the column from data to plot, usually 'C_bc_rem_rel' or 'k_bc_reld'
    - exclude_obs: can be used to exclude specific observations (e.g. the ones corresponding to non-pyrolysed biomass)
    - titles: list of str, titles for each subplot
    - others matplotlib parameters

    '''
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    if len(obs_sets)>1:
        axl = axes.ravel()
    else:
        axl = [axes]
    
    if isinstance(color_list, str):
        L = list(metadata[color_list].unique())
        L.remove('nan')
        color_dict = dict(zip(L, sns.color_palette("hls", 6)))

        # Initialize empty lists for legend handles and labels
        legend_handles = []
        legend_labels = []

    for i, s in enumerate(obs_sets):
        for j, o in enumerate(s):
            #print(j, o)
            if o in exclude_obs:
                continue
            sub_df = select_mineralization_data(ID_obs= o, data=data)
            if(len(sub_df) == 0):
                continue
            if(not sub_df[col_to_plot].isnull().values.any()):
                pivot = pd.pivot_table(sub_df, columns=['ID_obs'], values=[col_to_plot],  index=['time'])*factor
                pivot['time'] = pivot.index
                # select color based on criteria
                if isinstance(color_list, str):
                    # str is passed: e.g. BiomassClass, colorise then by biomass class
                    cat = metadata[color_list][o]
                    cc = np.array([color_dict[cat]])
                    # plot scatter
                    scatter = pivot.plot.scatter(ax=axl[i], x='time', y=col_to_plot,
                                    marker=marker_list[j // len(color_list)], c=cc, label=cat)
                    
                    # Add scatter plot to legend only if the category is not already in the legend
                    if cat not in legend_labels:
                        proxy = mpatches.Patch(color=color_dict[cat], label=cat)
                        legend_handles.append(proxy)
                        legend_labels.append(cat)

                    
                else:
                    # color_list is a palette, just random colors
                    cc = np.array([color_list[j % len(color_list)]])
                    # plot scatter
                    pivot.plot.scatter(ax=axl[i], x='time', y=col_to_plot,
                                    marker=marker_list[j // len(color_list)], c=cc)
                if line:
                    # plot a thin line; for decay rates, do not plot t=0, where decay = 0 for all obs
                    pivot.loc[1:,:].plot(kind='line', ax=axl[i], x='time', y=col_to_plot, color='grey', lw=lw, ls='dashed')
                    
            axl[i].set_xlabel('Time (days)', fontsize=14)
            axl[i].set_ylabel(ylabel, fontsize=14)
            axl[i].annotate(titles[i], xy=(0.02,1.02) , xycoords= 'axes fraction', color='black', fontsize=14)

            if isinstance(color_list, str):
                # Create color legend
                legend = axl[i].legend(handles=legend_handles, labels=legend_labels, title=legend_title, fontsize=16, title_fontsize=18)

    fig.subplots_adjust(top=0.93, wspace=0.2, hspace=0.27)
    if saveFig:
        plt.savefig(pathFig, dpi=dpiFig, bbox_inches='tight')

    return fig, axes, axl

## Curve fitting visualisation
def sort_params(f_model, params):
    '''
    Sort parameters from the dataframe, so they match the order of the params in the model function.
    
    Inputs:
    - f_model : a python function, with arguments (t, p1, p2, p3)
    - params : a dataframe

    Returns: 
    - f_params: list of param names (str) in the order of f_model 
    - p_opt_values: values of the params, in same order as required by f_model 
    '''
    param_names = params.index.get_level_values(1).to_list()
    param_values = params.values
    dic = {p:v for p,v in zip(param_names, param_values)}
    f_params = getfullargspec(f_model)[0][1:] # excluding time
    p_opt_values =[dic[k] for k in f_params]
    return f_params, p_opt_values

def select_multi_fits(ID_obs, fitdata, data, metadata, 
                      fit_limit=None, fit_rank='bic', fit_ascending=True, 
                      exclude_non_physical=False,
                      col_to_plot ='C_bc_rem_rel'):
    '''
    Prepare the data, for plotting multiple fits. This includes:
    - selecting the time series data for a given `ID_obs`, and `col_to_plot ='C_bc_rem_rel'`
    - performs physical checks on the fitted data, if `exclude_non_physical=True`
    - ranks fits according to `fit_rank`, and `fit_ascending=True`

    '''   
    # incubation data
    x, y = select_timeseries(ID_obs, data, col_to_plot, factor=100)
    xdataHighRes = np.arange(0,x[-1]+1)

    # fit data
    tpl_ID_obs = ('observation', 'ID_obs', 'int')
    tpl_rank = ('goodness', fit_rank, '1')
    tpl_bic = ('goodness', 'bic', '1')
    
    if isinstance(fitdata.columns, pd.MultiIndex):
        sub = fitdata[fitdata[tpl_ID_obs] == ID_obs]
        sub = sub.dropna(subset=[tpl_bic]) # remove failed fits, i.e. rows where all params are NaN or when rank is nan
        sub = sub.sort_values(tpl_rank, ascending=fit_ascending) # sort values by ascending bic
    else:
        sub = fitdata[fitdata['ID_obs'] == ID_obs]
        sub = sub.dropna(subset=['bic']) 
        sub = sub.sort_values(fit_rank, ascending=fit_ascending) 

    # exclude non-physical solutions, based on couple of checks
    if exclude_non_physical:
        col_check = [('checks', 'decay rates positive', 'bool'), 
                     ('checks', 'pool sizes positive', 'bool'),
                     ('checks', 'pool sizes below100', 'bool')]
        mask = sub[col_check]
        sub = sub[mask.all(axis=1)]

    # limit to x best fit, catching error if fit_limit > length 
    if fit_limit is not None:
        sub = sub.iloc[0:fit_limit, :] if fit_limit < len(sub.index) else sub.iloc[:, :] # select subset of fits      
        
    return ID_obs, x, y, xdataHighRes, sub

def plot_multi_fits(ID_obs, x, y, xdataHighRes, sub, metadata, ax=None, fig=None, saveTo=False):
    '''
    Plots a set fitting, for a given observation given by its ID_obs, alongside the decomposition data. Plot is added to a given axis, if passed.

    Parameters:
    - fit_limit: how many fit to plot among available; if None, plots all available, if less (int) e.g. 5 only plots 5 based on the fit_rank parameter
    - fit_rank: default `bic`, ranks the fits from lowest to highest `bic`, the lowest, the better (must be one of the goodness estimator (r2, chisqr, redchisqr, aic, bic)
    
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))
    
    ## prepare the data
    if isinstance(sub.columns, pd.MultiIndex):
        AuthorDate = metadata.loc[ID_obs, 'AuthorDate']
        BiomassClass = metadata.loc[ID_obs, 'BiomassClass']
        HHT=metadata.loc[ID_obs, 'HHT']
        HC= metadata.loc[ID_obs, 'H/C_tot'] if not (metadata.loc[ID_obs, 'H/C_tot'] == 'na' ) else -1
        tpl_model = ('strategy', 'model', 'str')
        tpl_method = ('strategy', 'method', 'str')
        tpl_r2 = ('goodness', 'r2', '1')
        tpl_ID_obs = ('observation', 'ID_obs', 'int')
        tpl_bic = ('goodness', 'bic', '1')
        tpl_DW = ('goodness', 'dw', '1')
    else:
        AuthorDate = metadata.loc[ID_obs, 'AuthorDate']
        BiomassClass = metadata.loc[ID_obs, 'BiomassClass']
        HHT=metadata.loc[ID_obs, 'HHT']
        HC= metadata.loc[ID_obs, 'H/C_tot'] if not (metadata.loc[ID_obs, 'H/C_tot'] == 'na' ) else -1
        tpl_model = 'model'
        tpl_method = 'method'
        tpl_r2 = 'r2'
        tpl_ID_obs = 'ID_obs'
        tpl_bic = 'bic'
        tpl_DW = 'dw'

    ## display the data
    # incubation data
    ax.plot(x, y, 'k.', ms=7, zorder=100,
        label='data for Obs={obs}, {AuthorDate}, {Bio}, HHT={HHT:.0f}, H/C={HC:.2f}'.format(
            obs=ID_obs, AuthorDate=AuthorDate, Bio=BiomassClass,
            HHT=HHT,
            HC= HC
        ))
       
    # fit data
    for i, fit in sub.iterrows():
        if ('strategy', 'strategy', 'str') in fit.index:
            strat = fit[('strategy', 'strategy', 'str')]
        else:
            strat = 'fit'
        f_model = map_model_function[fit[tpl_model]]
        method = fit[tpl_method]
        r2= fit[tpl_r2]
        bic=fit[tpl_bic]
        dw = fit[tpl_DW]
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
        p_opt_names, p_opt = sort_params(f_model, params)
        
        ax.plot(xdataHighRes, f_model(xdataHighRes, *p_opt), '-',
                label='{fit}, R\u00B2= {r2:.5f}, BIC= {bic:.2f}, DW= {dw:.2f}'.format(
                        fit=strat+': '+f_model.__name__+'_'+str(method), 
                        r2=r2, 
                        bic=bic, 
                        dw=dw),
                       )
    ## annotate the chart
    ax.set_xlabel('Time (days)', fontsize=16)
    ax.set_ylabel('Biochar C remaining (%)', fontsize=16) 
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(fontsize=14);
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax


def plot_multi_fits_decay(ID_obs, x, y, xdataHighRes, sub, metadata, factor=1, ax=None, fig=None, saveTo=False):
    '''
    Plots a set fitting, for a given observation given by its ID_obs, alongside the decomposition data. Plot is added to a given axis, if passed.

    Parameters:
    - fit_limit: how many fit to plot among available; if None, plots all available, if less (int) e.g. 5 only plots 5 based on the fit_rank parameter
    - fit_rank: default `bic`, ranks the fits from lowest to highest `bic`, the lowest, the better (must be one of the goodness estimator (r2, chisqr, redchisqr, aic, bic)
    
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))

    ## prepare the data
    if isinstance(sub.columns, pd.MultiIndex):
        AuthorDate = metadata.loc[ID_obs, 'AuthorDate']
        BiomassClass = metadata.loc[ID_obs, 'BiomassClass']
        HHT=metadata.loc[ID_obs, 'HHT']
        HC= metadata.loc[ID_obs, 'H/C_tot'] if not (metadata.loc[ID_obs, 'H/C_tot'] == 'na' ) else -1
        tpl_model = ('strategy', 'model', 'str')
        tpl_method = ('strategy', 'method', 'str')
        tpl_r2 = ('goodness', 'r2', '1')
        tpl_ID_obs = ('observation', 'ID_obs', 'int')
        tpl_bic = ('goodness', 'bic', '1')
        tpl_DW = ('goodness', 'dw', '1')
    else:
        AuthorDate = metadata.loc[ID_obs, 'AuthorDate']
        BiomassClass = metadata.loc[ID_obs, 'BiomassClass']
        HHT=metadata.loc[ID_obs, 'HHT']
        HC= metadata.loc[ID_obs, 'H/C_tot'] if not (metadata.loc[ID_obs, 'H/C_tot'] == 'na' ) else -1
        tpl_model = 'model'
        tpl_method = 'method'
        tpl_r2 = 'r2'
        tpl_ID_obs = 'ID_obs'
        tpl_bic = 'bic'
        tpl_DW = 'dw'

    ## display the data
    # incubation data
    ax.plot(x[1:], y[1:], 'k.', ms=7, zorder=100,
        label='data for Obs={obs}, {AuthorDate}, {Bio}, HHT={HHT:.0f}, H/C={HC:.2f}'.format(
            obs=ID_obs, AuthorDate=AuthorDate, Bio=BiomassClass,
            HHT=HHT,
            HC= HC
        ))
       
    # fit data
    for i, fit in sub.iterrows():
        f_model = map_model_function[fit[tpl_model]]
        method = fit[tpl_method]
        r2= fit[tpl_r2]
        bic=fit[tpl_bic]
        dw = fit[tpl_DW]
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
        p_opt_names, p_opt = sort_params(f_model, params)
        
        ax.plot(xdataHighRes[1:], -factor*np.diff(f_model(xdataHighRes, *p_opt)), '-',
                label='fit: {fit}, R\u00B2= {r2:.5f}, BIC= {bic:.2f}, DW= {dw:.2f}'.format(
                        fit=f_model.__name__+'_'+str(method), 
                        r2=r2, 
                        bic=bic, 
                        dw=dw),
                       )
    ## annotate the chart
    ax.set_yscale('log')
    ax.set_xlabel('Time (days)', fontsize=16)
    ax.set_ylabel(r"Biochar C decay rate (gC gC$_{0}^{-1}$ day$^{-1}$)", fontsize=16) 
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(fontsize=14);
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax
    

def plot_multi_fits_extended(ID_obs, x, y, xdataHighRes, sub, metadata, TH=101, ax=None, fig=None, saveTo=False):
    '''
    Plots the decay curves etrapolated to 100 years. Plot is added to a given axis, if passed.
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))
    
    ## prepare the data
    xdata100 = np.arange(0, 365*TH, 10)
    tpl_model = ('strategy', 'model', 'str')

    ## display the data
    # incubation data
    ax.plot(x, y, 'k.', ms=7, zorder=100)
    
    # fit data
    for i, fit in sub.iterrows():
        f_model = map_model_function[fit[tpl_model]]
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
        p_opt_names, p_opt = sort_params(f_model, params)
        ax.plot(xdata100, f_model(xdata100, *p_opt), '-')
    
    ## annotate the chart
    x_ticks = [day for day in xdata100 if day%3650 == 0] if TH < 500 else [day for day in xdata100 if day%36500 == 0] 
    x_labels = [str(i//365) for i in x_ticks]
    #ax2.yaxis.tick_right()
    ax.set_xlabel('Time (years)', fontsize=16)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    ax.set_xlim(0, 365*TH)
    ax.set_ylabel('Biochar C remaining (%)', fontsize=16)
    ax.set_ylim(0, 101)
    ax.tick_params(axis='both', which='major', labelsize=14)
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax

def plot_multi_fits_extended_decay(ID_obs, x, y, xdataHighRes, sub, metadata, TH=101, factor=1, ax=None, fig=None, saveTo=False):
    '''
    Plots the decay curves etrapolated to 100 years. Plot is added to a given axis, if passed.
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))
    
    ## prepare the data
    xdata100 = np.arange(0, 365*TH, 1)
    tpl_model = ('strategy', 'model', 'str')

    ## display the data
    # incubation data
    ax.plot(x[1:], y[1:], 'k.', ms=7, zorder=100)
    
    # fit data
    for i, fit in sub.iterrows():
        f_model = map_model_function[fit[tpl_model]]
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
        p_opt_names, p_opt = sort_params(f_model, params)
        ax.plot(xdata100[1:], -1*factor*np.diff(f_model(xdata100, *p_opt)), '-')
    
    ## annotate the chart
    x_ticks = [day for day in xdata100 if day%3650 == 0] if TH < 500 else [day for day in xdata100 if day%36500 == 0] 
    x_labels = [str(i//365) for i in x_ticks]
    #ax2.yaxis.tick_right()
    ax.set_yscale('log')
    ax.set_xlabel('Time (years)', fontsize=16)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    ax.set_xlim(0, 365*TH)
    ax.set_ylabel(r"Biochar C decay rate (gC gC$_{0}^{-1}$ day$^{-1}$)", fontsize=16) 
    ax.tick_params(axis='both', which='major', labelsize=14)
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax        

def plot_multi_fits_residuals(ID_obs, x, y, xdataHighRes, sub, metadata, ax=None, fig=None, saveTo=False):
    '''
    Plots the residuals of the fits. Plot is added to a given axis, if passed.
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))
    
    ## prepare the data
    tpl_model = ('strategy', 'model', 'str')
    tpl_residuals = ('goodness', 'residuals', 'vector')
    
    ## display the data
    # fit data
    for i, fit in sub.iterrows():
        f_model = map_model_function[fit[tpl_model]]
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
        p_opt_names, p_opt = sort_params(f_model, params)
        # get the residuals, saved as list without commas...
        t = fit[tpl_residuals]
        t = re.sub("\["+"\s+", "[", t)
        t = re.sub("\s+", ",", t)
        t = literal_eval(t)  
        ax.plot(f_model(x, *p_opt), t, '.', markersize=12)
    
    ## annotate the chart
    ax.axhline(y=0, xmin=0, xmax=1, color='black')
    ax.set_ylabel('Residuals', fontsize=14)
    ax.set_xlabel('Modelled biochar C remaining (%)', fontsize=14)

    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax  

def plot_multi_fits_report(ID_obs, x, y, xdataHighRes, sub, metadata, columns=None, ax=None, fig=None, saveTo=False):
    '''
    Plots the residuals of the fits. Plot is added to a given axis, if passed.
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(15,10))

    
    ## prepare the data
    if columns is None:
        # default set of columns to display
        columns = [# fit strategy
                     ('strategy', 'model', 'str'),
                     ('strategy', 'library', 'int'),
                     ('strategy', 'method', 'str'),
                     ('strategy', 'bounds', 'tuple'),
                     # goodness of fit
                     ('goodness', 'bic', '1'),
                     ('goodness', 'r2', '1'),
                     ('goodness', 'dw', '1'),
                     # estimator
                    ('estimator', 'BC100', '%'),
                    ('estimator', 'MRTa', '1'),
                    ('estimator', 't12', '1'),
                    # estimator q10c
                    ('q10correction', 'BC100c', '%'),
                    # checks
                #    ('checks', 'decay rates positive', 'bool'),
                #    ('checks', 'pool sizes positive', 'bool'),
                #    ('checks', 'pool sizes below100', 'bool'),
                #    ('checks', 'parameter constrained', 'bool'),
                    ]
    table = sub[columns].reset_index(drop=True)
    table = table.droplevel(level=[0,2], axis=1)
    table.update(table.loc[:, table.dtypes == float].applymap('{:,.2f}'.format))

    ## display the data
    ax.axis('off')
    font_size=16
    bbox=[0, 0, 1, 1]
    mpl_table = ax.table(cellText = table.values, rowLabels = table.index, bbox=bbox, colLabels=table.columns,)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)
    mpl_table.auto_set_column_width(col=list(range(len(table.columns)))) # Provide integer list of columns to adjust

    ## annotate the chart
    ax.set_title("Table", fontsize=16, weight='bold')
    fig.tight_layout()
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax  

def plot_multi_fits_snapshots(ID_obs, x, y, xdataHighRes, sub, metadata, snapshots=[100, 200, 500], fs=18, ax=None, fig=None, saveTo=False):
    '''
    Plots amount of biochar C remaining after times given as years in snapshots. Plot is added to a given axis, if passed.
        
    e.g. snapshots=[100, 200, 500] > 100 years, 200 years, 500 years
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))

    ## prepare the data
    tis = [k*365 for k in snapshots]
    tpl_model = ('strategy', 'model', 'str')
    Ts = metadata.loc[ID_obs, 'IncubationTemperature']
    
    ## display the data
    for i, fit in sub.iterrows():
        f_model = map_model_function[fit[tpl_model]]
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
        p_opt_names, p_opt = sort_params(f_model, params)
        bcs = [f_model(ti, *p_opt) for ti in tis]
        ax.scatter(x=snapshots, y=bcs, marker='X')
    
    ## annotate the chart
    ax.set_xlabel('Time (years)', fontsize=16)
    ax.set_ylabel('Biochar C remaining (%)', fontsize=16)
    ax.set_ylim(0, 101)
    ax.tick_params(axis='both', which='major', labelsize=14)
    note = "Soil at "+str(Ts)+"\u00B0C"
    ax.annotate(note, xy=(0.35, 0.90), xycoords="axes fraction", fontsize=fs)
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax  

def plot_multi_fits_snapshots_q10(ID_obs, x, y, xdataHighRes, sub, metadata, snapshots=[100, 200, 500], tTs= 14.9, fs=18, ax=None, fig=None, saveTo=False):
    '''
    Plots amount of biochar C remaining after times given as years in snapshots. Plot is added to a given axis, if passed.

        e.g. snapshots=[100, 200, 500] > 100 years, 200 years, 500 years
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))

    ## prepare the data
    tis = [k*365 for k in snapshots]
    tpl_model = ('strategy', 'model', 'str')
    Ts = metadata.loc[ID_obs, 'IncubationTemperature']
    
    # calculate Q10, fT for specific Ts pair
    def Q10(Ts, tTs=tTs):
        '''calculate Q10 factor, based on relationship given Woolf 2021 ES&T
        TO_DO: alternative: provide Q10 data, and calculate a relationship
        '''
        return 1.1 + ( 63.1579*(np.exp(-0.19*tTs) - np.exp(-0.19*Ts))/(Ts-tTs) )
    
    def fT(Ts, tTs=tTs):
        '''calculate the fT ratio of the decay rates at exp temperature Ts over target temperature tTs, based on Woolf 2021 ES&T'''
        return np.exp( np.log(Q10(Ts, tTs)) * (tTs-Ts)/10 )
    
    if Ts == tTs:
        # Nothing to change - plot same as before
        fig, ax = plot_multi_fits_snapshots(ID_obs, x, y, xdataHighRes, sub, metadata, snapshots, fs, ax, fig, saveTo)
    
    else:
        q10 = Q10(Ts, tTs)
        ft = fT(Ts, tTs)
        
        ## display the data
        for i, fit in sub.iterrows():
            f_model = map_model_function[fit[tpl_model]]
            params = fit.iloc[sub.columns.get_level_values(0)=='parameter'].dropna()
            p_opt_names, p_opt = sort_params(f_model, params)
            p_opt_c = []
            for i, p in enumerate(p_opt):
                if 'k' in p_opt_names[i]:
                    p_opt_c.append(p*ft)
                else:
                    p_opt_c.append(p)
            if fit[tpl_model] == 'powerModel':
                bcs = [powerModel_q10(ti, ft, *p_opt_c) for ti in tis]
            else:
                bcs = [f_model(ti, *p_opt_c) for ti in tis]

            ax.scatter(x=snapshots, y=bcs, marker='X')

        ## annotate the chart
        ax.set_xlabel('Time (years)', fontsize=16)
        ax.set_ylabel('Biochar C remaining (%)', fontsize=16)
        ax.set_ylim(0, 101)
        ax.tick_params(axis='both', which='major', labelsize=14)
        note = "Soil at "+str(tTs)+"\u00B0C \n(Q10 correction)"
        ax.annotate(note, xy=(0.35, 0.70), xycoords="axes fraction", fontsize=fs)

        # save chart
        if saveTo is not False:
            print("saving chart to:", saveTo)

    return fig, ax  

def rebuild_cov(f_model, params, stdev, covar):
    '''
    Rebuild the covariance matrix, in correct order from data saved as Excel...
    '''    
    f_params = getfullargspec(f_model)[0][1:] # excluding time

    param_names = params.index.get_level_values(1).to_list()
    param_values = params.values
    param_std_values = stdev.values # in same order as above
    covar = covar.droplevel(level=[0,2])
    
    def search_param(param_names, p):
        '''Returns an integer, position of search of param in params'''
        for s in range(len(params)):
            if p in param_names[s]:
                return s
    
    def search_covar(covar, p1, p2):
        '''Returns an integer, position of search covariance in matrix'''
        for s in range(len(covar.index)):
            if p1 in covar.index[s] and p2 in covar.index[s]:
                return s

    dicVal = {p:v for p,v in zip(param_names, param_values) if p in f_params}
    dicStdev = {p:v for p,v in zip(param_names, param_std_values) if p in f_params}
    #print(param_names, param_values, param_std_values)
    #print("dicStdev \n", dicStdev)
    
    p_opt_values =[dicVal[k] for k in f_params]
    
    p_cov = np.zeros((len(f_params),len(f_params)))
    for i in range(len(f_params)):
        for j in range(len(f_params)):
            if i == j:
                p_cov[i,j] = dicStdev[f_params[i]]**2
            else:
                p1 = f_params[i]
                p2 = f_params[j]
                s = search_covar(covar, p1, p2)
                p_cov[i,j] = covar[s]
    
    return f_params, p_opt_values, p_cov

def plot_best_fits_uncertainty(ID_obs, x, y, xdataHighRes, sub, metadata, bic_limit=3, ax=None, fig=None, saveTo=False):
    '''
    Plots the selected best fits, with their uncertainty ranges. Plot is added to a given axis, if passed.
    
    - bic_limit=3 >> only plots the 3 best fits, ranked by lowest bic 
    
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))

    
    ## prepare the data
    xdata100 = np.arange(0, 365*101, 10)
    tpl_model = ('strategy', 'model', 'str')

    ## display the data
    # fit data
    for i, fit in sub.iloc[:bic_limit, :].iterrows():
        f_model = map_model_function[fit[tpl_model]]
        
        params = fit.iloc[sub.columns.get_level_values(0)=='parameter']##.dropna() # ok to drop na, as we should not have na params
        stdev = fit.iloc[sub.columns.get_level_values(0)=='stddev_abs']##.dropna() # drops too many... here only.. cases where stdev could not be estimated are lost
        covar = fit.iloc[sub.columns.get_level_values(0)=='covariance']##.dropna() # drops too many... here only.. cases where covariance could not be estimated are lost
        ##p_opt_names, p_opt = sort_params(f_model, params)
        p_opt_names, p_opt, p_cov = rebuild_cov(f_model, params, stdev, covar)
        
        ax.plot(xdata100, f_model(xdata100, *p_opt), '-')
        if len(plt.gca().lines)>0:
            clr = plt.gca().lines[-1].get_color() # color of the last trace added, i.e. line above
        else:
            clr='b'
            
        ax.plot(xdata100, 
               f_model(xdata100, *p_opt) + 1.96* std_fmodel(xdata100, p_opt, np.transpose(p_cov), f_model),
               '-.', color=clr)
        ax.plot(xdata100, 
               f_model(xdata100, *p_opt) - 1.96* std_fmodel(xdata100, p_opt, p_cov, f_model),
               '-.', color=clr)    
    
    
    ## annotate the chart
    x_ticks = [day for day in xdata100 if day%3650 == 0] # Only pull out full years
    x_labels = [str(i//365) for i in x_ticks]
    #ax2.yaxis.tick_right()
    ax.set_xlabel('Time (years)', fontsize=16)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    ax.set_xlim(0, 365*101)
    ax.set_ylabel('Biochar C remaining (%)', fontsize=16)
    #ax.set_ylim(0, 101)
    ax.tick_params(axis='both', which='major', labelsize=14)

    note = "Estimated uncertainty for each fit"
    ax.annotate(note, xy=(0.2, 0.95), xycoords="axes fraction", fontsize=16)

    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax  



def plot_multi_fits_template(ID_obs, x, y, xdataHighRes, sub, metadata, ax=None, fig=None, saveTo=False):
    '''
    Just a template function, serving as a base structure for plotting results from fits.
    
    Plot is added to a given axis, if passed.
    '''
    if ax == None:
        # create a local figure
        fig, ax = plt.subplots(figsize=(20,10))

    ## prepare the data
    
    ## display the data
    
    ## annotate the chart
    
    # save chart
    if saveTo is not False:
        print("saving chart to:", saveTo)
        
    return fig, ax  

## Comparison with Woolf 2021
def plot_woolf_comparison(ID_obs, subset_fitdata_ext, data, metadata, saveTo='simulations/woolf-comparison-p0/'):
    ID_obs, x, y, xdataHighRes, sub = select_multi_fits(ID_obs, subset_fitdata_ext, data, metadata)

    AuthorDate = metadata.loc[ID_obs, 'AuthorDate']
    BiomassClass = metadata.loc[ID_obs, 'BiomassClass']
    HHT=metadata.loc[ID_obs, 'HHT']
    HC= metadata.loc[ID_obs, 'H/C_tot'] if not (metadata.loc[ID_obs, 'H/C_tot'] == 'na' ) else -1

    #fig, ax = plt.subplots(figsize=(14,6))
    fig = plt.figure(figsize=(14,6))
    r = 1
    c = 3
    ax = plt.subplot(r,c,(1, 2)) # decay fits
    ax2 = plt.subplot(r,c,(3, 3)) # 100 year extrapo

    ax.plot(x, y, 'kX', ms=10, 
        label='data for Obs={obs}, {AuthorDate}, {Bio}, HHT={HHT:.0f}, H/C={HC:.2f}'.format(
            obs=ID_obs, AuthorDate=AuthorDate, Bio=BiomassClass,
            HHT=HHT,
            HC= HC
        ))

    map_model_function = {
        "singleExp":singleExp,
        "singleExp_u":singleExp_u,
        "doubleExp":doubleExp,
        "doubleExp_u":doubleExp_u,
        "tripleExp":tripleExp,
        "tripleExp_u":tripleExp_u,
        "powerModel":powerModel,
    }
    tpl_model = 'model'
    tpl_method = 'method'
    tpl_r2 = 'r2'
    tpl_ID_obs = 'ID_obs'
    tpl_bic = 'bic'
    tpl_DW = 'dw'

    sub = sub.iloc[0] # there is only one row here! 
    f_model = map_model_function[sub[tpl_model]]
    method = sub[tpl_method]
    r2= sub[tpl_r2]
    bic=sub[tpl_bic]
    dw = sub[tpl_DW]
    param_cols = ['k', 'c', 'k1', 'k2', 'c1', 'c2', 'k3', 'c0', 'b', 'm', 'c3',]
    params = sub[param_cols].dropna()

    def sort_params(f_model, params):
        '''
        Sort parameters from the dataframe, so they match the order of the params in the model function
        f_model : a python function, with arguments (t, p1, p2, p3)
        param_names : a list of str p1, p2, p3 
        param_values: a list of values, in the order of the param_names
        Returns: param_names and param_values in same order as required by f_model 
        '''
        param_names = params.index.to_list()
        param_values = params.values #.flatten()
        dic = {p:v for p,v in zip(param_names, param_values)}
        f_params = getfullargspec(f_model)[0][1:] # excluding time
        p_opt_values =[dic[k] for k in f_params]
        return f_params, p_opt_values

    p_opt_names, p_opt = sort_params(f_model, params)

    ax.plot(xdataHighRes, f_model(xdataHighRes, *p_opt), '-',
            label='fit: {fit}, R\u00B2= {r2:.5f}, BIC= {bic:.2f}, DW= {dw:.2f}'.format(
                    fit=f_model.__name__+'_'+str(method), 
                    r2=r2, 
                    bic=bic, 
                    dw=dw),
                   )
    
    xdata100 = np.arange(0, 365*101, 10)
    ax2.plot(xdata100, f_model(xdata100, *p_opt), '-',
            label='fit: {fit}'.format(fit=f_model.__name__+'_'+str(method)))
    # Woolf's fit to plot

    
    param_woolf_pool = ['C1','C2', 'C3']
    param_woolf_decay = ['k1_W', 'k2_W', 'k3_W']
    param_woolf = param_woolf_pool + param_woolf_decay
    params = sub[param_woolf].dropna()
    params[param_woolf_pool] = params[param_woolf_pool].multiply(other=100)
    params[param_woolf_decay] = params[param_woolf_decay].multiply(other=1/365)
    params.rename({'C1':'c1','C2':'c2', 'C3':'c3', 'k1_W':'k1', 'k2_W':'k2', 'k3_W':'k3'}, inplace=True)

    p_opt_names, p_opt = sort_params(tripleExp_u, params)

    if len(params[params == 0]) == 2:
        m = 'doubleExp_u'
    elif len(params[params == 0]) == 0:
        m = 'tripleExp_u'
    elif len(params[params == 0]) == 4:
        m = 'singleExp'
    ax.plot(xdataHighRes, tripleExp_u(xdataHighRes, *p_opt), '-',
            label='fit: {} from Woolf2021'.format(m)
            )
    ax2.plot(xdata100, tripleExp_u(xdata100, *p_opt), '-',
            label='fit: {} from Woolf2021'.format(m)
            )

    ## annotate the chart
    ax.set_xlabel('Time (days)', fontsize=16)
    ax.set_ylabel('Biochar C remaining (%)', fontsize=16) 
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(fontsize=14);

    x_ticks = [day for day in xdata100 if day%3650 == 0] # Only pull out full years
    x_labels = [str(i//365) for i in x_ticks]
    #ax2.yaxis.tick_right()
    ax2.set_xlabel('Time (years)', fontsize=16)
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_labels)
    ax2.set_xlim(0, 365*101)
    ax2.set_ylabel('Biochar C remaining (%)', fontsize=16)
    ax2.set_ylim(0-5, 100+5)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    if saveTo is not False:
        #print("saving chart to:", saveTo)
        fig.tight_layout()
        fig.savefig(saveTo+'ID_obs'+str(ID_obs)+'.png', dpi=150)        
        plt.close()
    return fig, ax

## PCA - biplot

def plot_pca_biplot(projection, components, features, c1, c2, limfeat=5, Y=None, fig=None, ax=None, cbartxt=False, annotate=False, texts=None):
    '''
    For 2 given components c1 and c2 (numbered from 0), plots a PCA scatter plot & a loading plot (arrows, with max limfeat sorted by longer arrows).
    
    Y = 
    
    More info:
    https://ostwalprasad.github.io/machine-learning/PCA-using-python.html
    https://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pca.aspx
    '''
    projection = projection[:,[c1,c2]] # shape 110 x 2, i.e. observations x selected components >> used for scatter
    components = np.transpose(components[[c1,c2],:]) # shape 33 x 2, i.e. features x selected components >> used for loading
    
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(10,8))
    
    # SCATTER
    xs = projection[:,0] # scatter
    ys = projection[:,1] # scatter   
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    xss = xs * scalex
    yss = ys * scaley
    im = ax.scatter(xss, yss, c=Y, edgecolor='none', alpha=0.5, cmap=plt.cm.get_cmap('jet', 10))
    if cbartxt is not False:
        cbar = fig.colorbar(mappable=im, ax=ax)
        cbar.set_label(cbartxt, rotation=270, labelpad=20)
    if annotate:
        for i, txt in enumerate(texts):
            ax.annotate(txt, xy=(xss[i], yss[i]), xycoords="data", fontsize=10)
            
    # LOADING
    n = components.shape[0] # nb features, 33
    df = pd.DataFrame(components)
    df['length'] = np.sqrt(df[0]**2+df[1]**2)
    df['features'] = features
    df = df.sort_values(by='length', axis=0, ascending=False)
    components = df.values
    features = df['features'].values
    for i in range(min(n, limfeat)):
        ax.arrow(0, 0, components[i,0], components[i,1],color = 'r',alpha = 0.5, head_width=0.01, head_length=0.02)
        if features is None:
            ax.text(components[i,0]* 1.15, components[i,1] * 1.15, "Var"+str(i+1), color = 'green', ha = 'center', va = 'center')
        else:
            ax.text(components[i,0]* 1.15, components[i,1] * 1.15, features[i], color = 'g', ha = 'center', va = 'center')
     
    ax.set_xlabel("PC"+str(c1+1))
    ax.set_ylabel("PC"+str(c2+1))

    return fig, ax


def plot_RF_model(rf, xs, ycols=0, ylabel='', fig=None, ax=None, figsize=(10,3), fX=None, fY=None):
    '''
    Plots the outcome of a trained random forest regressor (rf) over a given feature space (xs), selecting which target variable to plot with (ys).

    - If xs has 1 dimension (1 feature), then it's a simple 2D plot.
    - If xs has 2 dimensions, and both are continuous variables; then it's a 2D plot with permanence as color
    - If xs has 3 dimensions, and 1 of them is a discrete variable; then, we can plot as a 2D plot, for each value of the discrete variable
    ...
    
    Input:
    - rf: a sklearn regressor, with the rf.predict() method
    - xs: a list of np.arrays with values for each features, order of features must be same as in rf.features_names_in_ , 
    - ys: which target variables to select, if several ones
    - ylabel: list of corresponding text labels for the different target variables
    - fX, fY: the original datataset, for plotting as scatter onto the visualisation
    '''
    feature_names_in = rf.feature_names_in_ #list of feature names input
    
    if not isinstance(ycols, list):
        ycols = [ycols]
        ylabel = [ylabel]
    
    if isinstance(xs, np.ndarray):        
        if ax is None:
            fig, ax = plt.subplots(1,1, figsize=figsize)    
        ys = rf.predict(xs)
        [ax.plot(xs, ys[:,y], label=yl) for y, yl  in zip(ycols, ylabel)]
        ax.set_xlabel('feature: '+feature_names_in[0])
        ax.set_ylabel('prediction')
        plt.legend();
        if fX is not None:
            ax.scatter(fX, fY, alpha=0.5)
        
    if isinstance(xs, list) and len(xs) == 2:
        if ax is None:
            fig, ax = plt.subplots(1, len(ycols), figsize=figsize) # one plot for each target variable
        for i, (y, yl) in enumerate(zip(ycols, ylabel)):
            # heatmap and contours
            X,Y = np.meshgrid(xs[0], xs[1])
            Z = np.zeros(shape=X.shape)
            for j in range(X.shape[0]):
                mesh = pd.DataFrame({feature_names_in[0]:X[j], feature_names_in[1]:Y[j]})
                Z[j,:] = rf.predict(mesh)[:,y] # get only one column # shape 1-dim array, of length same as mesh.shape[0]
            
            im = ax[i].pcolormesh(X,Y,Z, cmap= plt.cm.get_cmap('jet', 10), alpha=0.7, vmin=0, vmax=100) #PiYG_r seismic_r
            cbar = fig.colorbar(mappable=im, ax=ax[i])
            cbar.set_label('prediction: '+yl, rotation=270, labelpad=20)
            ax[i].set_xlabel('feature: '+feature_names_in[0])
            ax[i].set_ylabel('feature: '+feature_names_in[1])

            if fX is not None:
                ax[i].scatter(fX.iloc[:,0] , fX.iloc[:,1], c=fY[i], alpha=0.9, cmap=plt.cm.get_cmap('jet', 10), vmin=0, vmax=100, 
                             edgecolors="white")
            #contour
            #class nf(float):
            #    def __repr__(self):
            #        s = f'{self:.1f}'
            #        return f'{self:.0f}' if s[-1] == '0' else s
            #CS = ax[i].contour(Y, X, Z, [0, 50, 80], colors=['red','red', 'red']) # Y, X, Z 
            #CS.levels = [nf(val) for val in CS.levels]

            # Label levels with specially formatted floats
            #if plt.rcParams["text.usetex"]:
            #    fmt = r'%r '
            #else:
            #    fmt = '%r'
            #ax[i].clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=14)
        fig
        fig.tight_layout(pad=1.4)
    
    if isinstance(xs, list) and len(xs) == 3:
        # if 1 is not a float...
        oneCat = False
        for c, x in enumerate(xs):
            if isinstance(x[0], str):
                oneCat = True
                cx = c
        if oneCat:
            if ax is None:
                fig, ax = plt.subplots(len(set(xs[cx])), len(ycols), figsize=figsize) # rows: each category, column: each target variable
            for u, c in enumerate(set(xs[cx])):
                for i, (y, yl) in enumerate(zip(ycols, ylabel)):
                    # heatmap and contours
                    xs.pop(cx)
                    X,Y = np.meshgrid(xs[0], xs[1])
                    Z = np.zeros(shape=X.shape)
                    for j in range(X.shape[0]):
                        mesh = pd.DataFrame({feature_names_in[0]:X[j], feature_names_in[1]:Y[j], feature_names_in[2]: np.repeat(c, X.shape[1])})
                        Z[j,:] = rf.predict(mesh)[:,y] # get only one column # shape 1-dim array, of length same as mesh.shape[0]
                    
                    ax[u, i].pcolormesh(X,Y,Z, cmap= plt.cm.get_cmap('jet', 10), alpha=0.7, vmin=0, vmax=100) #PiYG_r seismic_r
                    cbar = fig.colorbar(mappable=im, ax=ax[i])
                    cbar.set_label('prediction: '+yl, rotation=270, labelpad=20)
                    ax[u, i].set_xlabel('feature: '+feature_names_in[0])
                    ax[u, i].set_ylabel('feature: '+feature_names_in[1])

                    if fX is not None:
                        ax[u, i].scatter(fX.iloc[:,0] , fX.iloc[:,1], c=fY[i], alpha=0.9, cmap=plt.cm.get_cmap('jet', 10), vmin=0, vmax=100, 
                                     edgecolors="white")


        # if all are continuous
        print("Not yet developped")
    # for 3D - slices through 3rd dimension - volume plot
    # https://matplotlib.org/stable/gallery/event_handling/image_slices_viewer.html