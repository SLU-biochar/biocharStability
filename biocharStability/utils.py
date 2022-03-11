"""
-*- coding: utf-8 -*-

biochar stability / utils.py

set of utility functions to load, handle, save the data
"""


# imports
from logging import raiseExceptions
from operator import index
from setuptools_scm import meta

import itertools
from tabulate import tabulate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

import plotly.graph_objects as go # originally version 4.6; then 5.1, now updated to 5.4
import plotly.express as px
import seaborn as sns
import ternary as ter
from adjustText import adjust_text


# load metadata, data, inventory, schema
def load_metadata(filepath, sheet_name='metadata', index_col=None, usecols='A:BN', skiprows=3, schema_rows=6, last_row=None,):
    '''
    From the biochar stability database (.xlsx), loads the metadata sheets and parse the appropriate formats to the various columns.
    
    Returns:
    - pandas dataframe, metadata
    - pandas dataframe, schema
    - dictionary, metadata_column_sets (column names of metadata, grouped by categories (id, bio, pyr, incub, char) and numerical/non-numerical data type)

    Usage: \n
        metadata, schema, metadata_column_sets = load_metadata(filepath='biochar_incubation_database_2022-03-22_v0.xlsx', sheet_name='metadata', index_col=None, usecols='A:BN', skiprows=3, schema_rows=6, last_row=None)
    
    '''
    # datatypes of each column (metadata.dtypes)
    dts = {'ID_obs':int, 'ID_art':int,
        'AuthorDate':str, 'NameObs':str, 'Replicates':int, 'Biomass':str,
        'BiomassClass':str, 'BiomassLignin':float, 'Pyrolysis':str, 'PyrolysisClass':str,
        'HHT':float,'RT':float, 'HR':float, 'Soil type':str, 
        'Soil clay content':float, 'Soil sand content':float,'Soil silt content':float, 'Soil organic matter content':float,
        'SoilMoisture_Absolute':float, 'SoilMoisture_WHC':float, 'pH_soil':float, 'pH_bc-soil':float,
        'SoilOrigin':str, 'SoilGeolocation':str, 'LabField':str, 'Cultivated':str,
        'IncubationTemperature':float, 'IncubationDuration':float, 'ApplicationRate':float,
        'TotalFluxMeasurement':str, 'BiocharFluxDetermination':str,
        'OtherIncubationDescription':str, 'AssociatedControl':str,
        'Ash550':float, 'Ash700':float,
        'Carbon':float, 'Carbon, organic':float, 'Hydrogen':float, 'Nitrogen':float, 'Sulphur':float,
        'Oxygen':float, 'Magnesium':float, 'Potassium':float, 'Phosphorous':float, 'H/C_org':float, 'H/C_tot':float,
        'O/C_org':float, 'C/N':float, 'pH_H2O':float, 'pH_CaCl2':float, 'CEC':float, 'SA_N2':float, 'SA_CO2':float,
        'FixedCarbon':float, 'VolatileMatter':float, 'ParticleSize':float, 'BulkDensity':float,
        'TrueDensity':float, 'OxidationH2O2':float, 'BiocharYield':float,
        'δ ¹³C biochar':float, 'δ ¹³C soil':float,
        'ID_Lehmann2019':str, 'ID_Woolf2021':str, 'ID_Lehmann2021':str   
    }

    # subset of columns, numerical & str/categorical 
    identification = {'numerical':['ID_obs', 'ID_art'], 'other':['AuthorDate','NameObs','ID_Lehmann2019', 'ID_Woolf2021', 'ID_Lehmann2021' ]}
    biomass = {'numerical':['BiomassLignin', ], 'other':['Biomass', 'BiomassClass', ]}
    pyrolysis = {'numerical':['HHT', 'HR', 'RT', ], 'other':['Pyrolysis', 'PyrolysisClass', ]}
    biochar = {'numerical':['Ash550', 'Ash700', 'Carbon', 'Carbon, organic', 'Hydrogen', 'Nitrogen', 'Sulphur', 'Oxygen', 'Magnesium', 'Potassium', 'Phosphorous', 'H/C_org', 'H/C_tot', 'O/C_org', 'C/N', 'pH_H2O', 'pH_CaCl2', 'CEC', 'SA_N2', 'SA_CO2', 'FixedCarbon', 'VolatileMatter', 'ParticleSize', 'BulkDensity', 'TrueDensity', 'OxidationH2O2', 'BiocharYield' ],
            'other':[]}
    incubation = {'numerical':['Replicates', 'Soil clay content', 'Soil sand content', 'Soil silt content',
                            'Soil organic matter content', 'SoilMoisture_Absolute', 'SoilMoisture_WHC',
                            'pH_soil', 'pH_bc-soil', 'IncubationTemperature', 'IncubationDuration', 'ApplicationRate', 'δ ¹³C biochar', 'δ ¹³C soil'],
                'other':['Soil type', 'SoilOrigin', 'SoilGeolocation',
                        'LabField', 'Cultivated', 'TotalFluxMeasurement', 'BiocharFluxDetermination', 'AssociatedControl',  ]}
    all_numerical = biomass['numerical']+pyrolysis['numerical'] + biochar['numerical'] + incubation['numerical']
    
    metadata_column_sets = {'identification':identification, 'biomass':biomass, 'pyrolysis':pyrolysis, 'biochar':biochar, 'incubation':incubation, 'all_numerical':all_numerical}


    df = pd.read_excel(filepath, sheet_name=sheet_name, index_col=index_col, usecols=usecols, skiprows=skiprows)

    schema = df.iloc[0:schema_rows, :].copy(deep=True)
    schema.set_index('name', drop=True, append=False, inplace=True)
    
    if last_row is None:
        last_row = len(df)
    metadata = df.iloc[schema_rows:last_row, 1:].copy(deep=True)

    metadata = metadata.astype(dts, errors='ignore') 
        # ignore errors, as float nan cannot be converted to integer (some NaN lefts for papers where didn't collect data yet - non isotopic ones)
        # note: empty values are converted to NaN (type float, np.isnan) but manually inserted 'na' (not available) are kept as str
    metadata.set_index('ID_obs', drop=True, append=False, inplace=True) # set ID_obs as id, without resorting
    metadata.name = 'metadata'
    print("Metadata loaded, with {} rows".format(last_row))
    return metadata, schema, metadata_column_sets


def load_data(filepath, sheet_name='data', index_col=None, usecols='A:AA', skiprows=3, schema_rows=4, last_row=None,):
    '''
    From the biochar stability database (.xlsx), loads the data sheet
    
    Returns:
    - pandas dataframe, data

    Usage: \n
        data = load_metadata(filepath='biochar_incubation_database_2022-03-22_v0.xlsx')
    
    '''
    df = pd.read_excel(filepath, sheet_name=sheet_name, index_col=index_col, usecols=usecols, skiprows=skiprows)
    if last_row is None:
        last_row = len(df)
    data = df.iloc[schema_rows:last_row, 1:].copy(deep=True) 
    data.name = 'data'
    print("Data loaded, with {} rows".format(last_row))
    return data

def load_articles(filepath, sheet_name='articles', index_col=None, usecols='A:P', skiprows=3, schema_rows=3, last_row=None,):
    '''
    From the biochar stability database (.xlsx), loads the article sheet
    
    Returns:
    - pandas dataframe, articles

    Usage: \n
        articles = load_articles(filepath='biochar_incubation_database_2022-03-22_v0.xlsx')
    '''
    df = pd.read_excel(filepath, sheet_name=sheet_name, index_col=index_col, usecols=usecols, skiprows=skiprows)

    if last_row is None:
        last_row = len(df)
    articles = df.iloc[schema_rows:last_row, 1:].copy(deep=True) 
    articles.name = 'articles'
    print("Articles loaded, with {} rows".format(last_row))
    return articles

def load_validation(filepath, sheet_name='validation', index_col=None, usecols='A:DD', skiprows=3, schema_rows=1, last_row=None,):
    '''
    From the biochar stability database (.xlsx), loads the article sheet
    
    Returns:
    - pandas dataframe, articles

    Usage: \n
        validation = load_validation(filepath='biochar_incubation_database_2022-03-22_v0.xlsx')
    '''
    df = pd.read_excel(filepath, sheet_name=sheet_name, index_col=index_col, usecols=usecols, skiprows=skiprows)

    if last_row is None:
        last_row = len(df)
    validation = df.iloc[schema_rows:last_row, 1:].copy(deep=True)
    validation.name = 'validation'
    print("Validation loaded, with {} rows".format(last_row))

    return validation

def load_q10data(filepath='biocharStability/database/biochar_q10_dataset_2022-05-27.xlsx', 
                 sheet_name='q10', index_col=None, usecols='B:Q', skiprows=1, schema_rows=0, last_row=None,):
    '''
    From the biochar q10 database (.xlsx), loads the data sheet
    
    Returns:
    - pandas dataframe, data

    Usage: \n
        data = load_metadata(filepath='biochar_incubation_database_2022-03-22_v0.xlsx')
    
    '''
    # datatypes of each column (metadata.dtypes)
    dts = {'q10_DOI':'str', 'q10_AuthorDate':'str', 'q10_ID_art':'int', 'T_1':'float', 'T_2':'float', 'T_avg':'float', 'Q10':'float',
       'BC_C':'float', 'BC_H':'float', 'H/C_mol':'float', 'BC-class-HT':'str', 'BC_class':'str', 'BC_HHT':'float',
       'Mineralisation over 1 year (%C ini)':'float', 'Q10 at 20C':'float'
    }
    df = pd.read_excel(filepath, sheet_name=sheet_name, index_col=index_col, usecols=usecols, skiprows=skiprows)
    if last_row is None:
        last_row = len(df)
    data = df.iloc[schema_rows:last_row, 0:].copy(deep=True) 
    data = data.astype(dts, errors='ignore') 
    data.name = 'q10'
    print("Q10 data loaded, with {} rows".format(last_row))
    return data
# database statistics
def print_database_completion_stats(articles=None,data=None, metadata=None, validation=None, schema=None):
    '''
    Displays some basic statistics about the loaded dataframes, if passed as argument
    '''

    ## Information available / Completeness of dataset
    def count_cells_info(df):
        '''
        Calculate some basic stats about the passed dataframe, such as:
        - number of data cells
        - number of cells that are not empty or NaN
        - number of cells that are "na", i.e. information not available for the observation
        - number of cells that are available information about the observation

        '''
        nb_cells = len(df)*len(df.columns)
        nb_info_na = df.isin(['na']).sum().sum() # (na = not available is also valuable information)
        nb_info_a = df.count().sum() - nb_info_na # nb of info available  (not empty, not NaN, minus not available) 
        nb_nan = df.isna().sum().sum()

        check_sum = nb_cells - nb_info_na - nb_info_a - nb_nan # should be equal to 0
        if check_sum != 0:
            raise Exception("Completion statistics do not add up to 100% for table: ", df.name)

        return {'Number of data cells': nb_cells,
                'Number of cells with \'available\' information': nb_info_a,
                'Number of cells with \'not available\' information':nb_info_na,
                'Number of empty cells (NaN)':nb_nan,
        }

    def print_count_cells_info(report):
        '''
        Takes dictionary report from `count_cells_info` and prints it in, with % values
        '''
        v0 = report['Number of data cells']
        for k,v in report.items():
            print(k, ':', v, '(', "{:.1f}".format(v/v0*100) ,'% )')

        #print("\n")
    
    print("** Completion statistics **")
    if(articles is not None):
        a_report = count_cells_info(articles)
        print("\t** table: articles **")
        print_count_cells_info(a_report)

    if(data is not None):
        d_report = count_cells_info(data)
        print("\t** table: data **")
        print_count_cells_info(d_report)

    if(metadata is not None):
        m_report = count_cells_info(metadata)
        print("\t** table: metadata **")
        print_count_cells_info(m_report)

    if(validation is not None):
        v_report = count_cells_info(validation)
        print("\t** table: validation **")
        print_count_cells_info(v_report)


    def count_article_obs_in_df(df):
        '''
        df must have columns ID_art and/or ID_obs, or as index
        '''
        c = 'ID_art'
        if c in df.columns:
            print('\tarticles:', len(df[c].unique())) # list(df[c].unique())
        
        c = 'ID_obs'
        if c in df.columns:
            print('\tobservations:', len(df[c].unique())) # list(df[c].unique())  
        elif c in df.index.name:
            print('\tobservations:', len(df.index.unique())) # list(df[c].unique())

    print("\n** Number of articles & observations from which `data` was extracted **")
    count_article_obs_in_df(data)
    print("\n** Number of articles & observations from which `metadata` was extracted **")
    count_article_obs_in_df(metadata)
    print("\n** Number of articles & observations for which `validation` data was recorded **")
    count_article_obs_in_df(validation)

# metadata statistics
def stats_by_class(df, c='BiomassClass', display=False, style='github'):
    r = pd.DataFrame()
    r['Counts'] = df[c].value_counts()
    r['%'] = 100* df[c].value_counts() / df[c].count()
    
    if display:
        print(tabulate(r, headers = 'keys', tablefmt = style, floatfmt=(".0f", ".0f", ".1f") ))
        print('\n')
    return r


def print_metadata_stats(metadata):
    '''
    Displays some statistics about the metadata available in the database, including:
    - table, biomass classes
    - table, pyrolysis classes

    '''
    print("Biomass types included:")
    biomassclass = stats_by_class(metadata, 'BiomassClass', True)
    
    print("Pyrolysis types included:")
    pyrolysisclass = stats_by_class(metadata, 'PyrolysisClass', True)

# manipulation of incubation data
def select_mineralization_data(ID_obs, data):
    '''
    Selects timeseries from data dataframe, for given observations. Does not raise warning if requested ID is not present in the data
    
    Returns: 
    - copy of subset of dataframe
    '''
    if not isinstance(ID_obs, list):
        ID_obs = [ID_obs]
    selected = data[(data['ID_obs'].isin(ID_obs))].copy(deep=True)

    return selected

def select_timeseries(ID_obs, data, col_to_plot, factor=1):
    '''
    Select a specific timeseries, for a single observations for the given column.
    
    Returns:
    -  numpy arrays of x (time) and y (col_to_plot), useful for plotting or curve_fitting
    
    Also applies checks if the data is available for that observation, and raises exception otherwise.
    '''
    sub_df = select_mineralization_data(ID_obs, data=data)
    if(len(sub_df) == 0):
        raise Exception("No decay data for this observation. ID_obs = {}".format(ID_obs))
    if(not sub_df[col_to_plot].isnull().values.any()):
        pivot = pd.pivot_table(sub_df, columns=['ID_obs'], values=[col_to_plot],  index=['time'])
    else:
        raise Exception("No {} data for this observation. ID_obs = {}".format(col_to_plot, str(ID_obs)))
        
    xdata = np.array(pivot.index.values).flatten()
    ydata = np.array(pivot[col_to_plot].values).flatten()*factor
    
    return xdata, ydata

def classify_observations_by_duration(metadata, months_cutoffs=[30, 20, 12, 0]):
    '''
    Creates a list of sets, each set containing the ID_obs of the observations being longer than x months.
    - [60, 48 , 24, 12, 10, 0] # 6 groups
    - [30, 15, 12, 10, 0] # 6 groups

    USAGE: \n
        months_cutoffs, sets_incubations, sets_labels = bs.classify_observations_by_duration(metadata)

    ''' 
    sets_incubations = []
    sets_labels = []
    for i, m in enumerate(months_cutoffs): 
        s = set(metadata[metadata['IncubationDuration'] > m*30.5]['IncubationDuration'].index)
        if i>0:
            j=0
            while(j<i):
                s = s - sets_incubations[j]
                j=j+1
        sets_incubations.append(s)
        sets_labels.append("Incubations longer than {} months ({} obs)".format(str(m), str(len(s))))
        print("Incubation longer than", m, "months ---> ", len(s), " observations")
        #print(s)

    return months_cutoffs, sets_incubations, sets_labels

## Function to filter observations by on criteria in the metadata or data
## To be used when building a stability model, correlating metadata to data
def merge_dfs(dfs, saveTo=False):
    '''Merge a list of dataframes `dfs`, with same columns, using pd.concat, and re-index it all'''
    df = pd.concat(dfs, axis=0, ignore_index=True)
    if saveTo is not False:
        df.to_excel(saveTo)
    return df

def intersection(lst1, lst2):
    '''Returns the intersection of two lists'''
    return list(set(lst1) & set(lst2))

def apply_has_decay_data(obs, metadata, data, withPrint=True):
    '''
    Exclude from the passed list of observations, observations which do not have decay data
    
    Function will:
        - print the number of observations exlcuded.
        - return the new list of observations
        - raise error if the new list is empty
    '''
    # checks if ID_obs is present in the data, but not if data in given column is specified...
    s = list(data[(data['ID_obs'].isin(obs))]['ID_obs'].unique())
    s = intersection(obs, s)
    
    if withPrint:
        print("Applying: HasDecayData. Excluded {} observations, reducing from {} to {} remaining observations.".format( len(obs)-len(s), len(obs), len(s) ) )
    
    if len(s) == 0:
        raise Exception("No observations left to analyse!")
    
    return s

def apply_exclude_by_duration(obs, metadata, excl_shorter_than=None, excl_longer_than=None, withPrint=True):
    '''
    Exclude from the passed list of observations, observations with durations shorter/longer than the parameters specified (expressed in days)
    
    Function will:
        - print the number of observations exlcuded.
        - return the new list of observations
        - raise error if the new list is empty
    '''
    
    if (excl_shorter_than is not None) & (excl_longer_than is None):
        s = metadata[ metadata['IncubationDuration'] >= excl_shorter_than ].index
    elif (excl_shorter_than is None) & (excl_longer_than is not None):
        s = metadata[ metadata['IncubationDuration'] <= excl_longer_than ].index
    else:
        s = metadata[ (metadata['IncubationDuration'] >= excl_shorter_than ) & (metadata['IncubationDuration'] <= excl_longer_than) ].index
    s = list(s)
    s = intersection(obs, s)
    
    if withPrint:
        print("Applying: IncubationDuration. Excluded {} observations, reducing from {} to {} remaining observations.".format( len(obs)-len(s), len(obs), len(s) ) )
    
    if len(s) == 0:
        raise Exception("No observations left to analyse!")
    
    return s

def apply_exclude_by_numeric(obs, metadata, var='HHT', excl_lower_than=None, excl_higher_than=None, withPrint=True):
    '''
    Exclude from the passed list of observations, observations with `var` lower/higher than the parameters specified (expressed in same unit as in metadata)
    
    Function will:
        - print the number of observations exlcuded.
        - return the new list of observations
        - raise error if the new list is empty
    '''
    
    if (excl_lower_than is not None) & (excl_higher_than is None):
        s = metadata[ metadata[var] >= excl_lower_than ].index
    elif (excl_lower_than is None) & (excl_higher_than is not None):
        s = metadata[ metadata[var] <= excl_higher_than ].index
    else:
        s = metadata[ (metadata[var] >= excl_lower_than ) & (metadata[var] <= excl_higher_than) ].index
    net_excl = len(metadata.index) - len(s) # number of exclusions for the tested criteria

    s = list(s)
    s = intersection(obs, s)
    
    if withPrint:
        #print("Applying: {}. Excluded {} observations, reducing from {} to {} remaining observations ({} matched the criteria).".format(var, len(obs)-len(s), len(obs), len(s), net_excl ) )
        print("Applying: {}. Excluded {} observations, reducing from {} to {} remaining observations.".format(var, len(obs)-len(s), len(obs), len(s) ) )
    
    if len(s) == 0:
        raise Exception("No observations left to analyse!")
    
    return s

def apply_exclude_exactly(obs, excl=[], withPrint=True):
    '''
    Exclude exactly the observations passed in the list excl.
    
    Function will:
        - print the number of observations exlcuded.
        - return the new list of observations
        - raise error if the new list is empty
    '''

    s =  list(set(obs) - set(excl))   
    if withPrint:
        print("Applying: RemoveExactly. Excluded {} observations, reducing from {} to {} remaining observations.".format(len(obs)-len(s), len(obs), len(s) ) )    
    if len(s) == 0:
        raise Exception("No observations left to analyse!")
    
    return s

def apply_intersect_best_fits(obs, df_best):
    '''
    Returns a dataframe of the best fits to use for the modelling, based on a list of observation passed. 
    '''
    for_model = df_best[ df_best['ID_obs'].isin(obs)]
    if len(for_model) != len(obs):
        dif = set(obs).difference(set(df_best['ID_obs']))
        print("Following observations are in `data` but have no `best-fit`: ", dif)
    print("Number of observations available for the correlation model: ", len(for_model))
    
    return for_model


## New correlation models
def obs_per_feat(X, plot=True):
    '''
    Given a metadata df with a subset of relevant features, counts number of nan, and plots nb of observations available without nan as a function of features
    '''
    Xna = X.isna().sum()
    X = X[Xna.sort_values().index]
    saved = []
    Xc = set()
    colums = []
    x = []
    y = []
    for i, c in enumerate(X.columns):
        Xc = Xc.union(set(X[X[c].isna() == True].index))   
        saved.append({
            'columns': colums.append(c),
            'n_obs_complete': len(X) - len(Xc),
            'n_obs_nan': len(Xc),
            'set_nan': Xc,
        })
        x.append(i+1)
        y.append(saved[i]['n_obs_complete'])
    
    if plot:
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.scatter(x,y, s=5)
        ax.plot(x,y, linewidth=0.5)
        ax.set_xlabel('Nb of features')
        ax.set_ylabel('Nb observations complete')
    return saved, x, y

def featcombin(X, n, m=0):
    '''
    Given a metadata df where columns are relevant features, creates a list containing lists of all combination of features, without repetititons
    '''
    sample_list = list(X.columns)
    list_combinations = list()
    for i in range(m, n + 1):
        list_combinations += list(itertools.combinations(sample_list, i))
    return list_combinations

def dry2daf(df, x, ash='Ash'):
    '''
    Convert quantity x expressed in dry-basis to a dry-ash-free basis, with ash content specified.
    New columns: named as `x_daf`
    How are NaNs handled: NaNs in one colum remain NaNs 
    '''
    dff = df.copy()
    dff.replace('na', np.NaN, inplace=True)
    if ash == 'Ash':
        # the mean of both ash content is taken, if any gaps 
        dff['Ash'] = dff[['Ash550', 'Ash700']].mean(axis=1)
    else:
        dff['Ash'] = dff[ash]
    
    dff[x+'_daf'] = dff[x]/(1-dff['Ash'])
    return dff