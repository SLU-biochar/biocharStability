"""
-*- coding: utf-8 -*-

biochar stability / analyse.py

set of utility functions to analyse the data - i.e. performing curve fitting on the data timeseries, apply Q10 temperature correction, calculate uncertainties on curve fits, build correlation models, ...
"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
import matplotlib.cm as cm

from inspect import getfullargspec, signature
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from statsmodels.stats.stattools import durbin_watson
from lmfit import Parameters, minimize, fit_report
from statsmodels.stats.stattools import durbin_watson

import seaborn as sns
import ternary as ter
from adjustText import adjust_text

from sklearn.compose import make_column_selector as selector
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline

from sklearn.preprocessing import StandardScaler, RobustScaler, LabelEncoder, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score, roc_curve, explained_variance_score
from sklearn.ensemble import RandomForestRegressor

from .utils import *

# The model function, f(x, …). It must take the independent variable as the first argument and the parameters to fit as separate remaining arguments.
model_names = {
    'singleExp':'S0',
    'doubleExp':'D0',
    'tripleExp':'T0',
    
    'singleExp_u':'S',
    'doubleExp_u':'D',
    'tripleExp_u':'T',
    
    'powerModel':'POW'
}
## MODEL FUNCTIONS
def singleExp(t, k):
    return 100*np.exp(-k*t)

def singleExp_u(t, k, c):
    return c*np.exp(-k*t)

def doubleExp(t, k1, k2, c1):
    return c1*np.exp(-k1*t)+ (100-c1)*np.exp(-k2*t)

def doubleExp_u(t, k1, k2, c1, c2):
    return c1*np.exp(-k1*t)+ c2*np.exp(-k2*t)

def tripleExp(t, k1, k2, k3, c1, c2):
    return c1*np.exp(-k1*t)+ c2*np.exp(-k2*t) + (100-c1-c2)*np.exp(-k3*t)

def tripleExp_u(t, k1, k2, k3, c1, c2, c3):
    return c1*np.exp(-k1*t)+ c2*np.exp(-k2*t) + c3*np.exp(-k3*t)

def powerModel(t, c0, b, m):
    # as defined in Zimmerman 2011
    if m > -1:
        z = c0 - c0*np.exp(b)/(m+1)*(t**(m+1))
    if m == -1:
        z = c0*(1-np.exp(b)*np.log(t))
    if m < -1:
        z = c0*np.exp(b)/(m+1)*(t**(m+1)) + c0*np.exp(b)/(m+1)*(0.001**(m+1)) + c0
    return z

def powerModel_q10(t, fT, c0, b, m):
    # as defined in Zimmerman 2011, with own recalculation for a different temperature based on k --> k*fT + re-integrate
    #print(m)
    if m > -1:
        z = c0 - c0*fT*np.exp(b)/(m+1)*(t**(m+1))
    elif m == -1:
        z = c0*(1-np.exp(b)*fT*np.log(t))
    elif m < -1:
        z = c0*fT*np.exp(b)/(m+1)*(t**(m+1)) + c0*fT*np.exp(b)/(m+1)*(0.001**(m+1)) + c0
    else:
        z = np.nan
    return z 

def powerModel_lim(t, c0, b, m, L=10*365):
    '''
    This function cannot be used in curve fitting step. Instead, curve fitting should be done with bs.powerModel.
    Then, the fitted parameters on bs.PowerModel are used to calculate the extrapolation with the assumption that the decay rates decreases only for L years.
    
    - L = 10*365 # days
    
    '''
    if hasattr(t, "__len__"):
        # yes, t is an array of values
        t_before = np.array(np.where(t <= L)).flatten()
        t_after = np.array(np.where(t>L)).flatten()
        if m > -1:
            z_before = c0 - c0*np.exp(b)/(m+1)*(t_before**(m+1))
        if m == -1:
            z_before = c0*(1-np.exp(b)*np.log(t_before))
        if m < -1:
            z_before = c0*np.exp(b)/(m+1)*(t_before**(m+1)) + c0*np.exp(b)/(m+1)*(0.001**(m+1)) + c0

        kL = (z_before[-1] - z_before[-2])/1.
        zL = z_before[-1]
        
        z_after = zL + kL * (t_after - L)
        
        z = np.append(z_before, z_after) # append, not add 
        
    else:
        # no, t is a single value
        if t < L:
            if m > -1:
                z = c0 - c0*np.exp(b)/(m+1)*(t**(m+1))
            if m == -1:
                z = c0*(1-np.exp(b)*np.log(t))
            if m < -1:
                z = c0*np.exp(b)/(m+1)*(t**(m+1)) + c0*np.exp(b)/(m+1)*(0.001**(m+1)) + c0
        else:
            # t >= L 
            if m > -1:
                zL = c0 - c0*np.exp(b)/(m+1)*(L**(m+1))
                zLn = c0 - c0*np.exp(b)/(m+1)*((L-1)**(m+1))
            if m == -1:
                zL = c0*(1-np.exp(b)*np.log(L))
                zLn = c0*(1-np.exp(b)*np.log(L-1))
            if m < -1:
                zL = c0*np.exp(b)/(m+1)*(L**(m+1)) + c0*np.exp(b)/(m+1)*(0.001**(m+1)) + c0
                zLn = c0*np.exp(b)/(m+1)*((L-1)**(m+1)) + c0*np.exp(b)/(m+1)*(0.001**(m+1)) + c0
            
            kL = (zL-zLn)/1.
            
            z = zL + kL * (t - L)
        
    return z

def k_singleExp(t, k):
    return (100*(k)*np.exp(-k*t))/100

def k_singleExp_u(t, k, c):
    return (c*(k)*np.exp(-k*t))/100

def k_doubleExp(t, k1, k2, c1):
    return (c1*(k1)*np.exp(-k1*t)+ (100-c1)*(k2)*np.exp(-k2*t))/100

def k_doubleExp_u(t, k1, k2, c1, c2):
    return (c1*(k1)*np.exp(-k1*t)+ c2*(k2)*np.exp(-k2*t))/100

def k_tripleExp(t, k1, k2, k3, c1, c2):
    return ( c1*(k1)*np.exp(-k1*t)+ c2*(k2)*np.exp(-k2*t) + (100-c1-c2)*(k3)*np.exp(-k3*t) )/100

def k_tripleExp_u(t, k1, k2, k3, c1, c2, c3):
    return (c1*(k1)*np.exp(-k1*t)+ c2*(k2)*np.exp(-k2*t) + c3*(k3)*np.exp(-k3*t))/100

def k_powerModel(t, c0, b, m):
    # as defined in Zimmerman 2011
    return (-c0 * np.exp(b) * t**m)/100

map_model_function = {
    "singleExp":singleExp,
    "singleExp_u":singleExp_u,
    "doubleExp":doubleExp,
    "doubleExp_u":doubleExp_u,
    "tripleExp":tripleExp,
    "tripleExp_u":tripleExp_u,
    "powerModel":powerModel,
}

## MODEL FUNCTIONS FOR LMFIT -- NOT NEEDED ANYMORE
def singleExp_u_lmfit(params,x,y):
    k1 = params['k1']
    c1 = params['c1']
    y_fit = singleExp_u(x, k1, c1)
    return y_fit - y

def doubleExp_u_lmfit(params,x,y):
    k1 = params['k1']
    k2 = params['k2']
    c1 = params['c1']
    c2 = params['c2']
    #print(c1.value, c2.value, k1.value, k2.value)
    y_fit = doubleExp_u(x, k1, k2, c1, c2)
    #print(y_fit-y)
    return y_fit - y

# some curve fitting metrics & stats
def rsquare(f, xdata, ydata, popt):
    '''Computes R-square value for given fit

    Input: \n
        - f : model function, callable

        - xdata: array

        - ydata: array

        - popt: returned by curve_fit

    Returns: \n
        - r2

        - residuals
    '''
    residuals = ydata-f(xdata, *popt) # '*' needed to unpack it (kwargs)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    
    rs = 1- ss_res/ss_tot
    return rs, residuals

def lmfit_fitting_stats(f, xdata, ydata, popt):
    '''Computes the fitting statistics (chi-square, reduced chi square, AIC, BIC) for a given fit, as implemented in lmfit 
    (source: https://github.com/lmfit/lmfit-py/blob/master/lmfit/minimizer.py line 364)
    Input: \n
        - f : model function, callable

        - xdata: array

        - ydata: array

        - popt: returned by curve_fit

    Outputs: \n
        - chisqr

        - redchisqr

        - AIC

        - BIC
        
    Usage: \n
        residuals, chisqr, redchisqr, aic, bic = lmfit_fitting_stats(f, xdata, ydata, popt)
        
    '''
    residuals = ydata-f(xdata, *popt) # residuals
    chisqr = np.sum(residuals**2) # sum-squared residuals
    ndata = len(residuals) # nb of data points
    nvarys = len(signature(f).parameters)-1 # nb of param in f function
    nfree = ndata - nvarys # degree of freedom: ndata - nparameters in the f function (i.e. nb args - 1 (since 1st arg is time))
    redchisqr = chisqr / max(1, nfree)
    chisqr = max(chisqr, 1.e-250*ndata) # modify chisqr value, if chisqr is 0, to be safe in the log
    _neg2_log_likel = ndata * np.log(chisqr / ndata)
    aic = _neg2_log_likel + 2 * nvarys # AIC
    bic = _neg2_log_likel + np.log(ndata) * nvarys # BIC
    
    return residuals, chisqr, redchisqr, aic, bic
   
def pretty_covariance(f, p_cov):
    '''
    Save a pretty covariance matrix as a table, with header names, in wide (matrix) and long dataframes

    Usage:\n
        df_matrix, df_long = pretty_covariance(f, p_cov)
    '''
    params = list(signature(f).parameters)[1:] # list of parameters, popping out first element (time)
    df_matrix = pd.DataFrame(p_cov, index=params, columns=params)
    df_long = pd.DataFrame(df_matrix.stack(), columns=['covariance'])
    
    return df_matrix, df_long

def print_fitting_report(f, xdata, ydata, p_opt, p_cov):
    '''
    Similar to LMFIT report \n
        [[Fit Statistics]]
            # fitting method   = Nelder-Mead
            # function evals   = 491
            # data points      = 44
            # variables        = 4
            chi-square         = 3.87297297 
            reduced chi-square = 0.09682432 
            Akaike info crit   = -98.9273573
            Bayesian info crit = -91.7905988
        [[Variables]]
            k1:  1.1683e-05 +/- 5.1274e-07 (4.39%) (init = 0)
            k2:  4.47275357 +/- 1.84402220 (41.23%) (init = 0)
            c1:  97.6338853 +/- 0.07100370 (0.07%) (init = 0)
            c2:  2.36617997 +/- 0.31931764 (13.50%) (init = 0)
        [[Correlations]] (unreported correlations are < 0.100)
            C(k1, c1) = 0.729
            C(k2, c2) = 0.387
            C(c1, c2) = -0.222
            C(k1, c2) = -0.162
            C(k2, c1) = 0.110
    '''
    residuals, chisqr, redchisqr, aic, bic = lmfit_fitting_stats(f, xdata, ydata, p_opt)
    df_matrix, df_long = pretty_covariance(f, p_cov)
    
    print("[[Fit Statistics]]")
    print("    # data points: ", len(residuals))
    print("    # variables: ", len(p_opt))
    print("    chi-square: ", chisqr)
    print("    reduced chi-square: ", redchisqr)
    print("    Akaike info crit: ", aic)
    print("    Bayesian info crit: ", bic)
    print("[[Variables]]")
    # for loop - standard deviations, absolute values
    for i, p in enumerate(df_matrix.columns):
        print("    {}: {} +/- {} ({} %)".format(p, p_opt[i], np.sqrt(df_matrix.loc[p,p]), np.sqrt(df_matrix.loc[p,p])/p_opt[i]*100 ))
    print("[[Correlations]] (unreported correlations are < 0.100)")
    # for loop - correlations, normalised as in lmfit
    correls = {}
    for i, p1 in enumerate(df_matrix.columns):
        for p2 in df_matrix.columns[i+1:]: # double for loop, to avoid dupplicates
            if p1 != p2:
                corr = df_matrix.loc[p1,p2] / ( np.sqrt(df_matrix.loc[p1,p1]) * np.sqrt(df_matrix.loc[p2,p2]) )
                correls[p1+", "+p2]= corr
    sort_correls = sorted(correls.items(), key=lambda it: abs(it[1]))
    sort_correls.reverse()
    for name, corr in sort_correls:
        if (corr >= 0.1) or (corr < -0.10):
            print("    C({}) = {}".format(name, corr) )


### SCIPY CURVE FITTING
def do_the_fit(f_model=doubleExp, xdata=[0], ydata=[0], p0=None, method=None, bounds=(-1*np.inf, np.inf), showPlot=False, full_output=False, maxfev=None):
    ''' Does the fitting for one observation, and returns the parameters of interest, using scipy.curve_fit

    Inputs: \n
        - f_model: choose one of the model functions to fit on, e.g. linear, exponential, singleExp, doubleExp, powerModel 
        - xdata: timesteps array, in days
        - ydata: decay array
        - p0: initial guess for optimal parameters, in curve_vit, otherwise set to 1
        - method: algorithm for fitting in curve_fit, {lm, trf, dogbox},
        - bounds = upper, and lower limits arrays for fitting parameters, in tuple format (lower_array, upper_array), each array with length of parameters or scalar
                Here, parameters are either: doubleExp: 0<k1,k2<c1<100
                                             singleExp: 0<k<c<100
                                             linear:    0< m < b < 100 (if function defined as -1*m*x + b)
                TO DO: pass bounds with the function choice?
        
    Outputs:  \n
        - fitted paramters, 
        - standard deviation of parameters, aboslute and relative
        - covariance matrix of fitted parameters
        - scaled correlations between parameters
        - stability estimates as BC100, apparent MRT, half-life 
        - r2 of fit (whether relevant or not)
        - other statisical indicator, chi, ... 
        - residuals of fit 

        - outcome of simple checks (passed failed)

    Usage:  \n
        p_opt, p_cov, p_std, r2, stab_dict = do_the_fit(f_model=doubleExp, xdata, ydata, p0=None, method=None, bounds=(-1*np.inf, np.inf))
    
    '''
    # Fitting
    try:
        fitting_output = {}
        if full_output: # for testing on single observations
            p_opt, p_cov, infodict, mesg, ier = curve_fit(f=f_model, xdata=xdata, ydata=ydata, 
                                 p0=p0, sigma=None, method=method, bounds=bounds, full_output=full_output, maxfev=maxfev)
                                 # sigma=ydata*0.05 ## to add noise/error on ydata
            return p_opt, p_cov, infodict, mesg, ier

        else:
            p_opt, p_cov,  = curve_fit(f=f_model, xdata=xdata, ydata=ydata, 
                                 p0=p0, sigma=None, method=method, bounds=bounds)
                                 # sigma=ydata*0.05 ## to add noise/error on ydata
        df_matrix, _ = pretty_covariance(f_model, p_cov)

        p_names = getfullargspec(f_model)[0] # p_names[0] = time, p_names[1:] = all other parameters, in same order as p_opt
        for j, p in enumerate(p_names[1:]):
            if 'k' in p:
                unit = '%/day'
            elif 'c' in p:
                unit = '%'
            else:
                unit='na'
            fitting_output[('parameter', p, unit)] = p_opt[j]
            fitting_output[('stddev_abs', 'std_'+p, 'abs')] = np.sqrt(df_matrix.loc[p,p])
            fitting_output[('stddev_rel', 'std_'+p, 'rel')] = np.sqrt(df_matrix.loc[p,p])/p_opt[j]*100 

        for i, p1 in enumerate(df_matrix.columns):
            for p2 in df_matrix.columns[i+1:]: # double for loop, to avoid dupplicates
                if p1 != p2:
                    fitting_output[('covariance', 'cov('+p1+', '+p2+')', '1')] = df_matrix.loc[p1,p2]
                    fitting_output[('correlation', 'corr('+p1+', '+p2+')', '%')] = df_matrix.loc[p1,p2] / ( np.sqrt(df_matrix.loc[p1,p1]) * np.sqrt(df_matrix.loc[p2,p2]) )
                    

        # Other parameters
        f_opt = f_model(xdata, *p_opt) # fitted curve
        p_std = np.sqrt(np.diag(p_cov)) # Standard deviation of each parameter  
        r2, residuals = rsquare(f_model, xdata, ydata, p_opt) # R2
        residuals, chisqr, redchisqr, aic, bic = lmfit_fitting_stats(f_model, xdata, ydata, p_opt)
        DW = durbin_watson(residuals)
        fitting_output[('goodness', 'r2', '1')] = r2
        fitting_output[('goodness', 'chisqr', '1')] = chisqr
        fitting_output[('goodness', 'redchisqr', '1')] = redchisqr
        fitting_output[('goodness', 'aic', '1')] = aic
        fitting_output[('goodness', 'bic', '1')] = bic
        fitting_output[('goodness', 'dw', '1')] = DW

        fitting_output[('goodness', 'residuals', 'vector')] = residuals # vector, likely not to save in an Excel table
        
        # Physical checks
        def check_decay_positivity(p_names, p_opt):
            '''Verifies that all decay rates have a positive value'''
            check = True
            for j, p in enumerate(p_names[1:]):
                if 'k' in p:
                    if p_opt[j] < 0:
                        check = False
            return check
        def check_pool_positivity(p_names, p_opt):
            '''Verifies that all pool sizes have a positive value'''
            check = True
            for j, p in enumerate(p_names[1:]):
                if 'c' in p:
                    if p_opt[j] < 0:
                        check = False
            return check

        def check_pool_below100(p_names, p_opt):
            '''Verifies that all pool sizes have a value below 100%'''
            check = True
            for j, p in enumerate(p_names[1:]):
                if 'c' in p:
                    if p_opt[j] > 100.0:
                        check = False
            return check

        def check_initial_carbon(f_model, p_opt):
            '''Calculates the total C at time = 0. Ideally should be as close as possible to 100%'''
            return f_model(0, *p_opt)

        def check_stdev_constraint(p_names, p_opt, df_matrix):
            '''Verifies that all parameters have a abs stdev smaller than absolute value of param
                i.e. np.abs(stdev) < np.aps(param)
            '''
            check = True
            for j, p in enumerate(p_names[1:]):
                    stdev = np.sqrt(df_matrix.loc[p,p])
                    if np.abs(stdev) > np.abs(p_opt[j]):
                        check = False
            return check
        def check_k_odm_difference(p_names, p_opt, odm=1):
            '''For exponential models only, with multiple pools, verifies that the decay rates have at least n order of magnitude difference when expressed in %/year
            '''
            k_list = []
            for j, p in enumerate(p_names[1:]):
                if 'k' in p:
                    k_list.append(p_opt[j]*365) # conversion to %/year, from %/day
            if len(k_list)>1:
                for i in range(len(k_list)):
                    for j in range(i + 1, len(k_list)):
                        if abs(np.log10(k_list[i]) - np.log10(k_list[j])) < 1:
                            return False
            return True

        fitting_output[('checks', 'decay rates positive', 'bool')] = check_decay_positivity(p_names, p_opt)
        fitting_output[('checks', 'pool sizes positive', 'bool')] = check_pool_positivity(p_names, p_opt)
        fitting_output[('checks', 'pool sizes below100', 'bool')] = check_pool_below100(p_names, p_opt)
        fitting_output[('checks', 'initial total carbon', '%')] = check_initial_carbon(f_model, p_opt)
        fitting_output[('checks', 'parameter constrained', 'bool')] = check_stdev_constraint(p_names, p_opt, df_matrix)
        fitting_output[('checks', 'decay rates 1 odm diff', 'bool')] = check_k_odm_difference(p_names, p_opt)

        # Stability values
        t100 = 100 * 365 # 100 years in days
        BC100 = f_model(t100, *p_opt) if f_model(t100, *p_opt) > 0 else 0 # remaning after 100 years

        def f12(t):
            return f_model(t, *p_opt) - 50.
        t12 = fsolve(f12, x0=50*365) # half-life time, 50% C remaining
        t12 = t12[0]
        
        def fMRT(t):
            return f_model(t, *p_opt) - np.exp(-1)
        MRT_app = fsolve(fMRT, x0=50*365) # apparent MRT, np.exp(-1) = 36.8% C remaining
        MRT_app = MRT_app[0]
        
        fitting_output[('estimator', 'BC100', '%')] = BC100
        fitting_output[('estimator', 't12', '1')] = t12/365
        fitting_output[('estimator', 'MRTa', '1')] = MRT_app/365

        stab_dict = {'BC100': BC100, 
                     't12': t12/365, # in years
                     'MRTa': MRT_app/365 # in years
                    }

        if showPlot:
            plt.scatter(xdata,ydata,c='black')
            plt.xlabel('Time (days)')
            plt.ylabel('y')
            plt.plot(xdata, f_model(xdata, *p_opt), c='red', ls='-', lw=2)
            plt.show()

            print_fitting_report(f=f_model, xdata=xdata, ydata=ydata, p_opt=p_opt, p_cov=p_cov)

        return p_opt, p_cov, p_std, r2, stab_dict, residuals, fitting_output
        
    except RuntimeError as e:
        # RuntimeError: Optimal parameters not found: The maximum number of function evaluations is exceeded.
        print(e)
        return [], [], [], [], [], [], {}
    
    except ValueError as e:
        print(e)
        return [], [], [], [], [], [], {}   

## LMFIT CURVE FITTING
def do_the_lmfit(f_model=doubleExp_u, 
               xdata=[0], ydata=[0],
               method='trf',
               bounds=(0, 100),
               showPlot=False,
               ):
    '''
    Does the fitting for one observation, and returns the parameters of interest, using lmfit.minimize

    Inputs: \n
        - f_model: choose one of the model functions to fit on, e.g. linear, exponential, singleExp, doubleExp, powerModel 
        - xdata: timesteps array, in days
        - ydata: decay array
        - bounds = upper, and lower limits arrays for fitting parameters, in tuple format (lower_array, upper_array), each array with length of parameters or scalar
                Here, parameters are either: doubleExp: 0<k1,k2,c1<100
                                             singleExp: 0<k<c<100
                                             linear:    0< m < b < 100 (if function defined as -1*m*x + b)
    
    TODO: implement use of initial guess, via p0 - list of initial guesses, in params.add(ar, value=, vary=True)
    '''
    fitting_output = {}
    pargs = getfullargspec(f_model).args[1:] # removing the first one, which is always time, list
    params = Parameters()
    for ar in pargs:
        params.add(ar, min=bounds[0], max=bounds[1])

    def f_model_lmfit(params, x, y):
        L = [params[k] for k in params]
        return f_model(x, *L) - y

    fitted_params = minimize(f_model_lmfit, params, args=(xdata, ydata,), method=method)
    p_opt = [fitted_params.params[p].value for p in fitted_params.var_names]
    r2, residuals = rsquare(f_model, xdata, ydata, p_opt)
    DW = durbin_watson(residuals)

    for j, p in enumerate(fitted_params.var_names):
        if 'k' in p:
            unit = '%/day'
        elif 'c' in p:
            unit = '%'
        else:
            unit='na'
        fitting_output[('parameter', p, unit)] = fitted_params.params[p].value
        fitting_output[('stddev_abs', 'std_'+p, 'abs')] = fitted_params.params[p].stderr
        fitting_output[('stddev_rel', 'std_'+p, 'rel')] = fitted_params.params[p].stderr/fitted_params.params[p].value*100 if fitted_params.params[p].stderr is not None else np.nan
    
    if hasattr(fitted_params, 'covar'):
        for i, p1 in enumerate(fitted_params.var_names):
            for j, p2 in enumerate(fitted_params.var_names[i+1:]): # double for loop, to avoid dupplicates
                if p1 != p2:
                    fitting_output[('covariance', 'cov('+p1+', '+p2+')', '1')] = fitted_params.covar[i,j]
                    fitting_output[('correlation', 'corr('+p1+', '+p2+')', '%')] = fitted_params.covar[i,j] / ( np.sqrt(fitted_params.covar[i,i]) * np.sqrt(fitted_params.covar[j,j]) )

    
   
    fitting_output[('goodness', 'r2', '1')] = r2
    fitting_output[('goodness', 'chisqr', '1')] = fitted_params.chisqr
    fitting_output[('goodness', 'redchisqr', '1')] = fitted_params.redchi
    fitting_output[('goodness', 'aic', '1')] = fitted_params.aic
    fitting_output[('goodness', 'bic', '1')] = fitted_params.bic
    fitting_output[('goodness', 'dw', '1')] = DW

    fitting_output[('goodness', 'residuals', 'vector')] = residuals # vector, likely not to save in an Excel table
    
    # Physical checks
    def check_decay_positivity(p_names, p_opt):
        '''Verifies that all decay rates have a positive value'''
        check = True
        for j, p in enumerate(p_names):
            if 'k' in p:
                if p_opt[j] < 0:
                    check = False
        return check
    def check_pool_positivity(p_names, p_opt):
        '''Verifies that all pool sizes have a positive value'''
        check = True
        for j, p in enumerate(p_names):
            if 'c' in p:
                if p_opt[j] < 0:
                    check = False
        return check

    def check_pool_below100(p_names, p_opt):
        '''Verifies that all pool sizes have a value below 100%'''
        check = True
        for j, p in enumerate(p_names):
            if 'c' in p:
                if p_opt[j] > 100.0:
                    check = False
        return check

    def check_initial_carbon(f_model, p_opt):
        '''Calculates the total C at time = 0. Ideally should be as close as possible to 100%'''
        return f_model(0, *p_opt)
    
    def check_stdev_constraint(p_names, p_opt, fitted_params):
        '''Verifies that all parameters have a abs stdev smaller than absolute value of param
            i.e. np.abs(stdev) < np.aps(param)
        '''
        check = True
        for j, p in enumerate(p_names):
            stdev = fitted_params.params[p].stderr
            if stdev is not None:
                if np.abs(stdev) > np.abs(p_opt[j]):
                        check = False
            else:
                    check = False
        return check
    
    def check_k_odm_difference(p_names, p_opt, odm=1):
            '''For exponential models only, with multiple pools, verifies that the decay rates have at least n order of magnitude difference when expressed in %/year
            '''
            k_list = []
            for j, p in enumerate(p_names[1:]):
                if 'k' in p:
                    k_list.append(p_opt[j]*365) # conversion to %/year, from %/day
            if len(k_list)>1:
                for i in range(len(k_list)):
                    for j in range(i + 1, len(k_list)):
                        if abs(np.log10(k_list[i]) - np.log10(k_list[j])) < 1:
                            return False
            return True
        
    fitting_output[('checks', 'decay rates positive', 'bool')] = check_decay_positivity(fitted_params.var_names, p_opt)
    fitting_output[('checks', 'pool sizes positive', 'bool')] = check_pool_positivity(fitted_params.var_names, p_opt)
    fitting_output[('checks', 'pool sizes below100', 'bool')] = check_pool_below100(fitted_params.var_names, p_opt)
    fitting_output[('checks', 'initial total carbon', '%')] = check_initial_carbon(f_model, p_opt)
    fitting_output[('checks', 'initial total carbon', '%')] = check_initial_carbon(f_model, p_opt)
    fitting_output[('checks', 'parameter constrained', 'bool')] = check_stdev_constraint(fitted_params.var_names, p_opt, fitted_params)
    fitting_output[('checks', 'decay rates 1 odm diff', 'bool')] = check_k_odm_difference(fitted_params.var_names, p_opt)

    # Stability values
    t100 = 100 * 365 # 100 years in days
    BC100 = f_model(t100, *p_opt) if f_model(t100, *p_opt) > 0 else 0 # remaning after 100 years

    def f12(t):
        return f_model(t, *p_opt) - 50.
    t12 = fsolve(f12, x0=50*365) # half-life time, 50% C remaining
    t12 = t12[0]
    
    def fMRT(t):
        return f_model(t, *p_opt) - np.exp(-1)
    MRT_app = fsolve(fMRT, x0=50*365) # apparent MRT, np.exp(-1) = 36.8% C remaining
    MRT_app = MRT_app[0]
    
    fitting_output[('estimator', 'BC100', '%')] = BC100
    fitting_output[('estimator', 't12', '1')] = t12/365
    fitting_output[('estimator', 'MRTa', '1')] = MRT_app/365

    # Pretty printing all the statistical data
    print(fit_report(fitted_params))
    
    L = [fitted_params.params[p].value for p in fitted_params.params] 
    if showPlot:
        plt.scatter(xdata,ydata,c='black')
        plt.xlabel('Time (days)')
        plt.ylabel('y')
        plt.plot(xdata, f_model(xdata, *L), c='red', ls='-', lw=2)
        plt.show()

    return fitting_output


def fit_all_observations(data, metadata, observations=[], fitting_strategies=[], library='scipy', variable='C_bc_rem_rel', factor=100, excel=''):
    '''
    Applies the fitting_strategies to all passed observations, using the scipy curve_fit function, and returns an extensive report as a list of dictionnary, optionalled saved as an Excel file

    USAGE: \n 
        df_scipy_fits = bs.fit_all_observations(data, metadata, fitting_strategies=fitting_strategies, library='scipy', variable='C_bc_rem_rel', factor=100, excel='scipy_allfits.xlsx') 
    '''
    all_fitting_outputs = []
    if len(observations) == 0:
        observations = metadata.index

    for j in observations:
        try:
            x, y = select_timeseries(j, data, variable, factor=factor)
        except Exception:
            continue
            
        for i, (f_model, method, p0, bounds) in enumerate(fitting_strategies):
            print(i, method, f_model, p0, bounds)
            if library == 'scipy':
                _, _, _, _, _, _, fitting_output = do_the_fit(f_model=f_model, xdata=x, ydata=y, p0=p0, method=method, bounds=bounds, showPlot=False)
            elif library == 'lmfit':
                fitting_output = do_the_lmfit(f_model=f_model, xdata=x, ydata=y, method=method, bounds=(0,1000), showPlot=False)
            else:
                raise Exception("Sorry, the fitting library is not know. Use either `scipy` or `lmfit` as keyworkd.")
            
            # adding fitting_strategy information to the fitting output
            fitting_output[('observation', 'ID_obs', 'int')] = j
            fitting_output[('strategy', 'id', 'int')] = i
            fitting_output[('strategy', 'library', 'int')] = library
            fitting_output[('strategy', 'method', 'str')] = method
            fitting_output[('strategy', 'model', 'str')] = f_model.__name__
            fitting_output[('strategy', 'bounds', 'tuple')] = bounds
            fitting_output[('strategy', 'p0', 'tuple')] = p0
            
            all_fitting_outputs.append(fitting_output)

    df = pd.DataFrame(all_fitting_outputs)
    L = df.columns
    new_col_order = []
    new_col_order.extend([ll for ll in L if ll[0] == 'observation'])
    new_col_order.extend([ll for ll in L if ll[0] == 'strategy'])
    new_col_order.extend([ll for ll in L if ll[0] == 'parameter'])
    new_col_order.extend([ll for ll in L if ll[0] == 'estimator'])
    new_col_order.extend([ll for ll in L if ll[0] == 'goodness'])
    new_col_order.extend([ll for ll in L if ll[0] == 'checks'])
    new_col_order.extend([ll for ll in L if ll[0] == 'stddev_abs'])
    new_col_order.extend([ll for ll in L if ll[0] == 'stddev_rel'])
    new_col_order.extend([ll for ll in L if ll[0] == 'correlation'])
    new_col_order.extend([ll for ll in L if ll[0] == 'covariance'])
    columns=pd.MultiIndex.from_tuples(new_col_order, names=['category', 'name', 'unit'])
    df = pd.DataFrame(df, columns=columns)
    df.index.name = 'index'
    if excel != '':
        df.to_excel(excel)

    return df

def load_fitted_observations(fp):
    '''
    Load as a dataframe, a previously saved set of fitting_outputs saved as an excel file
    '''
    return pd.read_excel(fp, index_col=0, header=[0,1,2])

## Q10 Soil Temperature Corrections
def applyQ10(fitdata, metadata, tTs=14.9, verbose=True):
    '''
    Applies the Q10 re-calibration of soil temperature to a calculated BC100 value at a given soil temperature to the 
    target soil temperature tTs
    
    Inputs: \n 
        - tTs: target soil temperature
        - metadata: dataframe with incubation meta-data
        - fitdata: dataframe containing all the fitting data, usually output of function `fit_all_observations`
    
    Outputs:  \n
        - a new dataframe based on fitdata, with a new column containing the recalibrated BC100 values.  
    '''
    
    def Q10(Ts, tTs=tTs):
        '''calculate Q10 factor, based on relationship given Woolf 2021 ES&T
        TO_DO: alternative: provide Q10 data, and calculate a relationship
        '''
        return 1.1 + ( 63.1579*(np.exp(-0.19*tTs) - np.exp(-0.19*Ts))/(Ts-tTs) )
    
    def fT(Ts, tTs=tTs):
        '''calculate the fT ratio of the decay rates at exp temperature Ts over target temperature tTs, based on Woolf 2021 ES&T'''
        return np.exp( np.log(Q10(Ts, tTs)) * (tTs-Ts)/10 )
    
    FIT_PARAMS = {
        'singleExp': {'ki': ['k'], 'ci':['']},
        'doubleExp': {'ki': ['k1', 'k2'], 'ci':['c1']},
        'tripleExp': {'ki': ['k1', 'k2', 'k3'], 'ci':['c1', 'c2']}, # k1, k2, k3, c1, c2
        
        'singleExp_u': {'ki': ['k'], 'ci':['c']},
        'doubleExp_u': {'ki': ['k1', 'k2'], 'ci':['c1', 'c2']},
        'tripleExp_u': {'ki': ['k1', 'k2', 'k3'], 'ci':['c1', 'c2', 'c3']},
        
        'powerModel': None, # Not supported correction factor (t, c0, b, m)'
    }
    
    df2 = fitdata.copy(deep=True)
    list_Q10 = []
    list_fT = []
    list_newBC100 = []
    list_oldBC100 = []
    
    for i, obs in fitdata.iterrows():
        ID_OBS = obs[('observation','ID_obs', 'int')]
        #print(i, ID_OBS)
        FIT_func = obs[('strategy','model', 'str')]
        #print(FIT_func)
        TS = metadata.loc[ID_OBS, 'IncubationTemperature']
        #print(TS)
        oldBC100 = obs[('estimator','BC100','%')]
        
        # find the fitting parameters in fitted_params
        fitparams = fitdata.xs('parameter', level='category', axis=1).loc[i]
        
        #cal Q10 and fT
        q10 = Q10(TS, tTs)
        list_Q10.append(q10)
        ft = fT(TS, tTs)
        list_fT.append(ft)
        
        # apply correction to decay constants
        if FIT_func == 'singleExp':
            newBC100 = singleExp(100*365, ft*fitparams['k'][0])
        elif FIT_func == 'singleExp_u':
            newBC100 = singleExp_u(100*365, ft*fitparams['k'][0], fitparams['c'][0])            
        elif FIT_func == 'doubleExp':
            newBC100 = doubleExp(100*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], fitparams['c1'][0])
        elif FIT_func == 'doubleExp_u':
            newBC100 = doubleExp_u(100*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], fitparams['c1'][0], fitparams['c2'][0])
        elif FIT_func == 'tripleExp':
            newBC100 = tripleExp(100*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], ft*fitparams['k3'][0], fitparams['c1'][0], fitparams['c2'][0])
        elif FIT_func == 'tripleExp_u':
            newBC100 = tripleExp_u(100*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], ft*fitparams['k3'][0], fitparams['c1'][0], fitparams['c2'][0], fitparams['c3'][0])
        elif FIT_func == 'powerModel':
            newBC100 = powerModel_q10(100*365, ft, fitparams['c0'][0], fitparams['b'][0], fitparams['m'][0], )
        else:
            # not supported
            newBC100 = obs[('estimator','BC100','%')] # copy previous value
            print("Warning: encoutered model function not supported by Q10 recalibration for obs: ", ID_OBS, FIT_func)
        list_newBC100.append(newBC100)
        list_oldBC100.append(oldBC100)
        #print(type(fitparams['k']))
        if verbose:
            print(i, len(fitparams), 'Q10: ', q10, 'fT: ', ft, FIT_func, oldBC100, newBC100)
        
    
        #if i > 2:
        #    break
    df2[('q10correction','Q10','1')] = list_Q10
    df2[('q10correction','fT','1')] = list_fT
    df2[('q10correction','BC100c','%')] = list_newBC100
    return df2 #, list_oldBC100, list_newBC100

def fit_Q10_data():
    '''
    Fits a relationship between Q10 and soil temperature, based on experimental data compiled.
    
    TODO: load Q10 data from excel, fit exponential model, return model parameters
    '''
    return None

## OTHER TEMPERATURE ADJUSTMENT METHODS (added 2023-07)

# Woolf2021
def Q10(T1, T2):
    '''calculate Q10 average factor, based on relationship given Woolf 2021 ES&T
    '''
    return 1.1 + ( 63.1579*(np.exp(-0.19*T1) - np.exp(-0.19*T2)) / (T2-T1) )

def fT_Woolf(T1, T2):
    '''calculate the fT ratio of the decay rates at exp temperature Ts over target temperature tTs, based on integration presented in Woolf 2021 ES&T'''
    return 1 if T1 == T2 else np.exp( np.log(Q10(T1, T2)) * (T2-T1)/10 )

# step-wise Woolf2021
def fT_WoolfStepwise(T1,T2, s=0.001):
    '''Calculates the fT with Woolf2021 method and a stepwise approach'''
    f=1
    sign = +1 if T1<T2 else -1
    x = np.arange(T1, T2, sign*s)
    for i in x:
        f=f*fT_Woolf(i,i+sign*s)
    return f

# Eye-fit on data
def fT_DataFit(T1, T2):
    '''from T1 to T2, using exponential relationship fitted on k_Y2 and Incubation temperature of whole dataset, by eye'''
    return (0.9 * np.exp(0.0200*T2) - 0.7) / (0.9 * np.exp(0.0200*T1) - 0.7)


def applyTempCorr(fitdata, metadata, tTs=14.9, TH=100, method=fT_Woolf):
    '''
    Applies the re-calibration of soil temperature to a calculated BC_TH value at a given soil temperature to the 
    target soil temperature tTs
    
    Inputs: \n 
        - tTs: target soil temperature, in degree C
        - TH: time horizon, in years
        - fT: pass a function used to calculate fT, among the fT availables:
                - fT_Woolf (or bs.fT_Woolf, if library is imported as bs)
                - fT_WoolfStepwise (or bs.)
                - fT_DataFit (or bs.)
        - metadata: dataframe with incubation meta-data
        - fitdata: dataframe containing all the fitting data, usually output of function `fit_all_observations`
    
    Outputs:  \n
        - a new dataframe based on fitdata, with a new column containing the recalibrated BC100 values.  
        - list_oldBC100: old values as list
        - list_newBC100: new values as list
    '''
    
    FIT_PARAMS = {
        'singleExp': {'ki': ['k'], 'ci':['']},
        'doubleExp': {'ki': ['k1', 'k2'], 'ci':['c1']},
        'tripleExp': {'ki': ['k1', 'k2', 'k3'], 'ci':['c1', 'c2']}, # k1, k2, k3, c1, c2
        
        'singleExp_u': {'ki': ['k'], 'ci':['c']},
        'doubleExp_u': {'ki': ['k1', 'k2'], 'ci':['c1', 'c2']},
        'tripleExp_u': {'ki': ['k1', 'k2', 'k3'], 'ci':['c1', 'c2', 'c3']},
        
        'powerModel': None, # Correction supported by separate function, powerModel_q10
    }
    
    df2 = fitdata.copy(deep=True)
    list_fT = []
    list_newBC100 = []
    list_oldBC100 = []
    
    for i, obs in fitdata.iterrows():
        ID_OBS = obs[('observation','ID_obs', 'int')]
        #print(i, ID_OBS)
        FIT_func = obs[('strategy','model', 'str')]
        #print(FIT_func)
        TS = metadata.loc[ID_OBS, 'IncubationTemperature']
        #print(TS)
        oldBC100 = obs[('estimator','BC100','%')]
        
        # find the fitting parameters in fitted_params
        fitparams = fitdata.xs('parameter', level='category', axis=1).loc[i]
        
        #calc fT
        ft = method(TS, tTs)
        list_fT.append(ft)
        
        # apply correction to decay constants
        if FIT_func == 'singleExp':
            newBC100 = singleExp(TH*365, ft*fitparams['k'][0])
        elif FIT_func == 'singleExp_u':
            newBC100 = singleExp_u(TH*365, ft*fitparams['k'][0], fitparams['c'][0])            
        elif FIT_func == 'doubleExp':
            newBC100 = doubleExp(TH*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], fitparams['c1'][0])
        elif FIT_func == 'doubleExp_u':
            newBC100 = doubleExp_u(TH*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], fitparams['c1'][0], fitparams['c2'][0])
        elif FIT_func == 'tripleExp':
            newBC100 = tripleExp(TH*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], ft*fitparams['k3'][0], fitparams['c1'][0], fitparams['c2'][0])
        elif FIT_func == 'tripleExp_u':
            newBC100 = tripleExp_u(TH*365, ft*fitparams['k1'][0], ft*fitparams['k2'][0], ft*fitparams['k3'][0], fitparams['c1'][0], fitparams['c2'][0], fitparams['c3'][0])
        elif FIT_func == 'powerModel':
            newBC100 = powerModel_q10(TH*365, ft, fitparams['c0'][0], fitparams['b'][0], fitparams['m'][0], )
        else:
            # not supported
            newBC100 = obs[('estimator','BC100','%')] # copy previous value
            print("Warning: encoutered model function not supported by Q10 recalibration for obs: ", ID_OBS, FIT_func)
        
        newBC100 = max(0, newBC100)
        list_newBC100.append(newBC100)
        list_oldBC100.append(oldBC100)

    #df2[('q10correction','Q10','1')] = list_Q10
    df2[('q10correction','fT','1')] = list_fT
    df2[('q10correction','BC_TH_'+str(TH)+'_TS_'+str(tTs),'%')] = list_newBC100
    return df2 , list_oldBC100, list_newBC100

## BEST FIT SELECTION
def select_best_fit(fitdata, model_pool=[], checks={}, spec_val={}, rank= ('bic', True), saveExcel=False):
    '''
    Based on the passed criteria, selects for all observations available fitted_data a single best fit if it exists. 

    The criterias are:
    1. `model_pool`: a list of model functions names (e.g. `singleExp`, `doubleExp`, `tripleExp`,...) that are considered, for the best fit (e.g. we can decide to exclude powerModel or singleExp)
    2. `checks`: a dict of ohysical checks that must be passed, to be considered a potential best fit, e.g. {'decay rates positive':TRUE}
    3. `spec_val`: any other key present in the fitdata data, e.g. {'library':'scipy'} to work only with scipy and exclude lmfit 
    4. `rank`: a tuple, among remaining fits, select the one with lowest/highest score in given goodness indicator, e.g. ('bic', True) means lowest BIC is best
    
    Inputs:
    - fitdata: a dataframe loaded from Excel, created by function `applyQ10` or `fit_all_observations` 
    - saveExcel: if False, does not save. Alternatively: give filepath to save to.
    
    Outputs:
    - dataframe, selected best fits for each observations, in same format as fitted_data; the dataframe is also saved as an Excel file (option)
    - a list of ID_obs where no best fit was found
    
    USAGE: \n
        no_best, df_best = select_best_fit(fitdata, model_pool, checks, spec_val)
        
    '''
    # fitted df - remove the sub-levels of multi-index, for easier access to properties
    fitdata = fitdata.droplevel(level=[0,2], axis=1)
    all_obs = set(fitdata['ID_obs'].unique())
    
    # apply model_pool filter
    if len(model_pool)>0:
        fitdata = fitdata[ fitdata['model'].isin(model_pool)]

    # apply checks filter
    if len(checks)>0:
        for k,v in checks.items():
            fitdata = fitdata[ fitdata[k] == v]
    
    # apply spec_val filter
    if len(spec_val)>0:
        for k,v in spec_val.items():
            fitdata = fitdata[ fitdata[k] == v]    
    
    # ranking and best fit selection 
    fitdata = fitdata.sort_values(by=['ID_obs', rank[0] ], ascending=[True, rank[1]] )
    df_best = fitdata.drop_duplicates('ID_obs')
    
    rem_obs = set(df_best['ID_obs'].unique())
    
    no_best = list(all_obs.difference(rem_obs))
    
    if saveExcel is not False:
        df_best.to_excel(saveExcel)
    
    return no_best, df_best

## ERROR PROPAGATION for specific models
def std_singleExp(t, p_opt, p_cov):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the model function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty
        
    See definition of `singleExp` for order of parameters
    '''
    k = p_opt[0]
    gX = np.array([-100*t*np.exp(-k*t) ]) # gradient of doubleExp_u, for vector X
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def std_singleExp_u(t, p_opt, p_cov):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the model function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    See definition of `singleExp_u` for order of parameters
    '''
    k = p_opt[0]
    c = p_opt[1]

    gX = np.array([-c*t*np.exp(-k*t) , np.exp(-k*t) ]) # gradient of doubleExp_u, for vector X
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)


def std_doubleExp(t, p_opt, p_cov):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the double exponential function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    See definition of `doubleExp_u` for order of parameters
    '''
    k1 = p_opt[0]
    k2 = p_opt[1]
    c1 = p_opt[2]
    gX = np.array([-c1*t*np.exp(-k1*t) , (100-c1)*t*np.exp(-k2*t) , np.exp(-k1*t) - np.exp(-k2*t) ]) # gradient of doubleExp_u, for vector X
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def std_doubleExp_u(t, p_opt, p_cov):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the double exponential function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    
    '''
    k1 = p_opt[0]
    k2 = p_opt[1]
    c1 = p_opt[2]
    c2 = p_opt[3]
    gX = np.array([-c1*t*np.exp(-k1*t) , - c2*t*np.exp(-k2*t) , np.exp(-k1*t) , np.exp(-k2*t) ]) # gradient of doubleExp_u, for vector X
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def std_tripleExp(t, p_opt, p_cov):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the triple exponential function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    
    '''
    k1 = p_opt[0]
    k2 = p_opt[1]
    k3 = p_opt[2]
    c1 = p_opt[3]
    c2 = p_opt[4]

    gX = np.array([-c1*t*np.exp(-k1*t) , - c2*t*np.exp(-k2*t) , (100-c1-c2)*t*np.exp(-k3*t) , np.exp(-k1*t) - np.exp(-k3*t)  , np.exp(-k2*t) - np.exp(-k3*t) ]) # gradient of doubleExp_u, for vector X
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def std_tripleExp_u(t, p_opt, p_cov):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the triple exponential function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    '''
    k1 = p_opt[0]
    k2 = p_opt[1]
    k3 = p_opt[2]
    c1 = p_opt[3]
    c2 = p_opt[4]
    c3 = p_opt[5]

    gX = np.array([-c1*t*np.exp(-k1*t) , - c2*t*np.exp(-k2*t) , - c3*t*np.exp(-k3*t) , np.exp(-k1*t) , np.exp(-k2*t) , np.exp(-k3*t) ]) # gradient of doubleExp_u, for vector X
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def std_powerModel(t, p_opt, p_cov, fT=1):
    '''
    Given optimal parameters of a fit (p_opt), the covariance matrix of the fitted parameters (p_cov),
    calculate the propagated error through the model function, based on Talor series expansion of 1st order
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    
    '''
    c0 = p_opt[0]
    b =  p_opt[1]
    m =  p_opt[2]
    if m > -1:
        gX = np.array([ 1 - fT * np.exp(b) / (m+1) * (t**(m+1)) , fT * c0 * np.exp(b) / (m+1) * (t**(m+1))  , -1 * fT * np.exp(b) / (m+1) / (m+1) * (t**(m+1)) * (np.log(t)*(m+1) - 1 ) ]) # gradient of doubleExp_u, for vector X
    if m == -1:
        gX = np.array([ 1 - fT * np.exp(b) * np.log(t) , - fT * c0 * np.log(t) * np.exp(b), 0]) # gradient of doubleExp_u, for vector X
    if m < -1:
        # taken same as for m > -1, since same equation but just a constant insert for time = 0 
        gX = np.array([ 1 - fT * np.exp(b) / (m+1) * (t**(m+1)) , fT * c0 * np.exp(b) / (m+1) * (t**(m+1))  , -1 * fT * np.exp(b) / (m+1) / (m+1) * (t**(m+1)) * (np.log(t)*(m+1) - 1 ) ]) # gradient of doubleExp_u, for vector X
    
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def std_fmodel(xdata, p_opt, p_cov, fmodel=doubleExp_u):
    '''
    Returns caculated propagated error for the given model, over the entire time range given (xdata)
    '''
    mapping = {
        singleExp:std_singleExp,
        singleExp_u:std_singleExp_u,
        doubleExp:std_doubleExp,
        doubleExp_u:std_doubleExp_u,
        tripleExp:std_tripleExp,
        tripleExp_u:std_tripleExp_u,
        powerModel:std_powerModel,
        } # add other functions here, when calculated (by hand)
    std_func = mapping[fmodel]  
    return np.array([std_func(x, p_opt, p_cov) for x in xdata])

## Linear Correlations - static
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, explained_variance_score

def rebuild_cov2(f_model, params, stdev, covar):
    '''
    Rebuild the covariance matrix, in correct order from data saved as Excel dataframe.
    '''    
    f_params = getfullargspec(f_model)[0][1:] # excluding time

    param_names = params.index.to_list()
    param_values = params.values
    param_std_values = stdev.values # in same order as above
    covar = covar #.droplevel(level=[0,2])
    
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

def apply_Q10_bis(Th, tTs, Ts, FIT_func, params, p_opt, p_cov):
    '''
    Given an observation, its model function, its parameters, a given TH and TS, the function:
    - calculates BC_TH_TS remaining at any given TH and TS
    - calculates a propagaged error on BC_TH_TS based on fitting uncertainty of fitted params
    
    IMPROVEMENT/TODO: at the moment, does not include uncertaity from the Q10 temperature correction (would need uncertainty in fit kT, Q10, as well as new derative calculations)
    
    '''
    
    def Q10(Ts, tTs=tTs):
        '''calculate Q10 factor, based on relationship given Woolf 2021 ES&T
        TO_DO: alternative: provide Q10 data, and calculate a relationship
        '''
        return 1.1 + ( 63.1579*(np.exp(-0.19*tTs) - np.exp(-0.19*Ts))/(Ts-tTs) )
    
    def fT(Ts, tTs=tTs):
        '''calculate the fT ratio of the decay rates at exp temperature Ts over target temperature tTs, based on Woolf 2021 ES&T'''
        return np.exp( np.log(Q10(Ts, tTs)) * (tTs-Ts)/10 )
    
    ft = 1 if Ts == tTs else fT(Ts, tTs)
    
    # apply correction to decay constants
    if FIT_func == 'singleExp':
        newBC100 = singleExp(Th*365, ft*params['k'])
        #sigmaBC100 = std_fmodel([Th*365*ft], [ft*params['k']], p_cov, singleExp)
        sigmaBC100 = std_fmodel([Th*365*ft], [params['k']], p_cov, singleExp)
    elif FIT_func == 'singleExp_u':
        newBC100 = singleExp_u(Th*365, ft*params['k'], params['c'])
        #sigmaBC100 = std_fmodel([Th*365*ft], [ft*params['k'], params['c']] , p_cov, singleExp_u)
        sigmaBC100 = std_fmodel([Th*365*ft], [params['k'], params['c']] , p_cov, singleExp_u)
    elif FIT_func == 'doubleExp':
        newBC100 = doubleExp(Th*365, ft*params['k1'], ft*params['k2'], params['c1'])
        #sigmaBC100 = std_fmodel([Th*365*ft], [ft*params['k1'], ft*params['k2'], params['c1']] , p_cov, doubleExp)
        sigmaBC100 = std_fmodel([Th*365*ft], [params['k1'], params['k2'], params['c1']] , p_cov, doubleExp)
    elif FIT_func == 'doubleExp_u':
        newBC100 = doubleExp_u(Th*365, ft*params['k1'], ft*params['k2'], params['c1'], params['c2'])
        #sigmaBC100 = std_fmodel([Th*365*ft], [ft*params['k1'], ft*params['k2'], params['c1'], params['c2']] , p_cov, doubleExp_u)
        sigmaBC100 = std_fmodel([Th*365*ft], [params['k1'], params['k2'], params['c1'], params['c2']] , p_cov, doubleExp_u)
    elif FIT_func == 'tripleExp':
        newBC100 = tripleExp(Th*365, ft*params['k1'], ft*params['k2'], ft*params['k3'], params['c1'], params['c2'])
        #sigmaBC100 = std_fmodel([Th*365*ft], [ft*params['k1'], ft*params['k2'], ft*params['k3'], params['c1'], params['c2']] , p_cov, tripleExp)
        sigmaBC100 = std_fmodel([Th*365*ft], [params['k1'], params['k2'], params['k3'], params['c1'], params['c2']] , p_cov, tripleExp)
    elif FIT_func == 'tripleExp_u':
        newBC100 = tripleExp_u(Th*365, ft*params['k1'], ft*params['k2'], ft*params['k3'], params['c1'], params['c2'], params['c3'])
        #sigmaBC100 = std_fmodel([Th*365*ft], [ft*params['k1'], ft*params['k2'], ft*params['k3'], params['c1'], params['c2'], params['c3']] , p_cov, tripleExp_u)
        sigmaBC100 = std_fmodel([Th*365*ft], [params['k1'], params['k2'], params['k3'], params['c1'], params['c2'], params['c3']] , p_cov, tripleExp_u)
    elif FIT_func == 'powerModel':
        newBC100 = powerModel_q10(Th*365, ft, params['c0'], params['b'], params['m'])
        sigmaBC100 = std_fmodel([Th*365*ft], [params['c0'], params['b'], params['m']] , p_cov, powerModel)
    else:
        # not supported
        print("Warning: encoutered model function not supported by Q10 recalibration for obs: ", FIT_func)
        newBC100 = np.nan
        sigmaBC100 = np.nan

    return newBC100, sigmaBC100

def correlation_linear(fitdata, metadata, x, y='BC_Th_Ts', Ts=15, Th=100, plot=True, interactive=False, ax=None, outliers=None, ifig=None, iobs=None, trendline=True, trenderror=False, wSigma=False, an=False, figsize=(9,6)):
    '''
    Performs a linear correlation between the two variables x and y (which must be one column in either fitdata or metadata),
    for a given soil temperature and a time horizon (if relevant, usually the case, when calculating BC_t). 
    
    Normal use case:
    - y='BC_Th_Ts' will calculate amount of biochar C remaining at Th and Ts passed as argument; while x is usually taken from the metadata columns
    - x must be from a metadata column, e.g. H/C_org, H/C_tot, HHT, Carbon, Carbon, organic, or another accepted value (currently H/C_all only)
        - If `x = H/C_all`, then the metadata considered will be: H/C_org whenever available, and H/C_tot to bridge gaps.

    - At time of plotting, the function will return the number of observations dropped due to missing data. 
    - Selected `outliers` can be removed by specifying a list of ID_obs 
    - Plotting can be done in static manner, or interactive mode (with sliders)
    
    '''
    fitdata = fitdata.copy(deep=True)
    obs = list(fitdata['ID_obs'].unique())
    # exclude identified outliers specified in *outliers*
    if outliers is not None:
        obs = set(obs)-set(outliers)

    fitdata = fitdata.set_index('ID_obs') # change index to ID_obs, so that series have the same index
    fitdata = fitdata[fitdata.index.isin(obs)].copy(deep=True) # filter metadata for the observations in the model minus outliers
    sub_metadata = metadata[metadata.index.isin(obs)] # filter metadata for the observations in the model
    sub_metadata = sub_metadata.replace('na',np.NaN) # replaces 'na' strings by NaN, for correct handling 
    
    # Make sure fitdata and sub_metadata have same index elements & order
    fitdata = fitdata.sort_index() # sort index, to be sure it matches with metadata
    sub_metadata = sub_metadata.sort_index() # sort index, to be sure it matches with metadata
    if not fitdata.index.equals(sub_metadata.index):
        # True if “other” is an Index and it has the same elements and order as the calling index; False otherwise.
        raise Exception("Warning: the indexes of fitdata & sub_metadata do not have same elments or same order")

    # Retrieving vX
    acc_x = ['H/C_all'] # special values of x that are accepted (e.g. H/C_all > using H/C_org if available, but bridging with H/C_tot if not available)
    if x not in metadata.columns and x not in acc_x: # verify x is in metadata.columns
        raise Exception("x = `{}` must be a column in `metadata` or an accepted value (check documentation)".format(x))
    elif x in metadata.columns:
        vX = sub_metadata[x]
    else: # x is a special value, swith and handle separately
        if x == 'H/C_all':
            vX = sub_metadata['H/C_org'].fillna(sub_metadata['H/C_tot'], inplace=False) 

    # Counting NaNs in vX
    vXnan = list(vX[vX.isna()].index)

    # Retrieving vY
    if y == 'BC_Th_Ts':
        # recalculate BC for given Th and Ts, and save it as vY, based on Q10 and fT saved in fitdata
        vY = []
        sigvY =[]
        for obs, row in fitdata.iterrows():
            Ts_ini = metadata.loc[obs, 'IncubationTemperature']
            FIT_func = row['model']
            params = row[['k', 'c', 'k1', 'k2', 'c1', 'c2', 'k3', 'c0', 'b', 'm', 'c3']]
            # std is in absolute value (other named .1 => relative value, in %)
            params_std = row[['std_k', 'std_c','std_k1', 'std_k2', 'std_c1', 'std_c2', 'std_k3', 'std_c3', 'std_c0', 'std_b', 'std_m']]
            params_corr = ['corr(k, c)', 'corr(k1, k2)', 'corr(k1, c1)', 'corr(k2, c1)', 'corr(k1, c2)', 'corr(k2, c2)', 'corr(c1, c2)', 'corr(k1, k3)',
                            'corr(k2, k3)', 'corr(k3, c1)', 'corr(k3, c2)', 'corr(k1, c3)', 'corr(k2, c3)', 'corr(k3, c3)', 'corr(c1, c3)', 'corr(c2, c3)','corr(c0, b)', 'corr(c0, m)', 'corr(b, m)']
            params_cov = row[['cov(k, c)', 'cov(k1, k2)', 'cov(k1, c1)', 'cov(k2, c1)', 'cov(k1, c2)', 'cov(k2, c2)',
                           'cov(c1, c2)', 'cov(k1, k3)', 'cov(k2, k3)', 'cov(k3, c1)','cov(k3, c2)', 'cov(k1, c3)', 'cov(k2, c3)', 'cov(k3, c3)','cov(c1, c3)', 'cov(c2, c3)', 'cov(c0, b)', 'cov(c0, m)', 'cov(b, m)']]
            map_model_function = {"singleExp":singleExp,"singleExp_u":singleExp_u,"doubleExp":doubleExp,"doubleExp_u":doubleExp_u,"tripleExp":tripleExp,"tripleExp_u":tripleExp_u,"powerModel":powerModel,}            
            
            _, p_opt, p_cov = rebuild_cov2(map_model_function[FIT_func], params, params_std, params_cov)
            BC, sBC = apply_Q10_bis(Th, Ts, Ts_ini, FIT_func, params, p_opt, p_cov)
            vY.append(BC)
            sigvY.append(sBC[0]) # [0] because functions returns a time series
            
            if obs == iobs:
                print("BC_TH_TS", BC)
                print("sigma_BC", sBC[0])
                
        fitdata['BC_Th_Ts'] = vY
        fitdata['std_BC_Th_Ts'] = sigvY
        vY = fitdata['BC_Th_Ts']
        sigvY = fitdata['std_BC_Th_Ts']
    elif y in fitdata.columns:
        vY = fitdata[y]
    else:
        raise Exception("y = `{}` must be a column in `fitdata` or an accepted value (check documentation)".format(y))

    # Counting NaNs in vY
    vYnan = list(vY[vY.isna()].index)

    # Remove mutual any observation with nan values
    vXYnan = set(vXnan + vYnan)
    vX = vX.drop(vXYnan, axis=0).sort_index(axis=0)
    vY = vY.drop(vXYnan, axis=0).sort_index(axis=0)
    sigvY = sigvY.drop(vXYnan, axis=0).sort_index(axis=0) # removing irrelevant observations, as in vX, vY
    sigvY = sigvY.fillna(0) # remainin nan, where uncertainty could not be calculated for some reason, replaced by 0

    # Names of Obs, for interactive plot; and colors
    names = ["Obs "+str(i) for i in vX.index]
    norm = plt.Normalize(1,4)
    cmap = plt.cm.RdYlGn

    vXa = np.array(vX).reshape(-1,1)
    vYa = np.array(vY).reshape(-1,1)
    sigvYa = 1.96*sigvY.values #np.array([sigvY.values, sigvY.values])
    xxx = np.array([np.min(vXa),np.max(vXa)])

    ## WITH SKLEARN
    # Linear Curve fitting - without testing/training sets - full descriptive analysis
    reg_hc = LinearRegression()
    reg_hc.fit(X = vXa, y = vYa )
    y_predict_reg_hc = reg_hc.predict(X = vXa )
    r2_regHC = reg_hc.score(X=vXa, y=vYa)
    rmse_regHC = np.sqrt(mean_squared_error(vYa, y_predict_reg_hc))
    exvar_regHC = explained_variance_score(vYa, y_predict_reg_hc)
    sl = float(reg_hc.coef_)
    it = float(reg_hc.intercept_)
    
    ## WITH SCIPY
    def linear(x, a, b):
        return a*x+b
    
    def std_linear(x, p_opt, p_cov):
        a = p_opt[0]
        b = p_opt[1]
        gX = np.array([x, 1])
        std_f2 = gX.T @ p_cov @ gX
        return np.sqrt(std_f2)
    
    p_opt1, p_cov1 = curve_fit(f=linear, xdata=vX, ydata=vY,p0=None, method='lm', sigma=None, absolute_sigma=False)
    r2_1, _ =  rsquare(linear, vX, vY, p_opt1) # R2
    if wSigma:
        p_opt2, p_cov2 = curve_fit(f=linear, xdata=vX, ydata=vY,p0=p_opt1, method='lm', sigma=1/sigvY, absolute_sigma=False)
        r2_2, _ =  rsquare(linear, vX, vY, p_opt2) # R2
        ## why 1/sigvY suggested by some (as weighing factor) and why not just sigvY (which does not find a solution?)

    df = pd.DataFrame(index=vX.index)
    df[x] = vX.values
    df[y] = vY.values
    df['BiomassClass'] = [v for o, v in sub_metadata['BiomassClass'].items() if o in vX.index]
    df['PyrolysisClass'] = [v for o, v in sub_metadata['PyrolysisClass'].items() if o in vX.index]

    # plotting
    if plot:
        if len(vXnan)>0:
            print("The parameter x=`{}` contains {} NaN values, from observations: {}".format(x, len(vXnan), vXnan))
        if len(vYnan)>0:
            print("The parameter y=`{}` contains {} NaN values, from observations: {}".format(y, len(vYnan), vYnan))
        if outliers is not None:
                print("Specified outliers ({}) were removed from regression: {}".format(len(outliers), outliers))

        print("Remaining number of observations in the graph:{}".format(len(vXa)))
        # print("R2 of H:C ratio regression: ", r2_regHC)
        # print("RMSE of H:C ratio regression: ", rmse_regHC)
        # print("Explained variance: ", exvar_regHC)

        name="Linear regression between: {} and {}".format(x, y)
        #fig, ax = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=[12,6])
        fig = plt.figure(figsize=figsize)
        r = 1
        c = 3
        ax1 = plt.subplot(r,c,(1, 3)) # decay fits
        #ax2 = plt.subplot(r,c,(3, 3)) # 100 year extrapo
        ax2=''
        ax = [ax1, ax2]
        i=0

        #ax[i].scatter(vXa, vYa, s=20, marker="o")
    
        sns.scatterplot(data=df, x=x, y=y,
                        style='PyrolysisClass', hue='BiomassClass', palette=sns.color_palette("Set2", as_cmap=False)[:len(df['BiomassClass'].unique())], 
                        ax = ax[i], s=80)
        ax[i].set_xlabel(x, fontsize=14)
        ax[i].set_ylabel(y, fontsize=14)
        ax[i].tick_params(axis='both', which='major', labelsize=12)
        ax[i].legend(facecolor='white', framealpha=0.2, fontsize=8) 

        #ax[i].errorbar(vXa.flatten(), vYa.flatten(), yerr=sigvYa, fmt="none",
        #               barsabove=True, elinewidth=0.5, ecolor='k', solid_capstyle='butt', capthick=0.5, capsize=2)

        if an is not False:
            texts = [str(ii) for ii in vX.index]
            xss = list(vXa)
            yss = list(vYa)
            for ii, txt in enumerate(texts):
                ax1.annotate(txt, xy=(xss[ii]+0.01, yss[ii]), xycoords="data", fontsize=6)

        #ax[i].plot(xxx, reg_hc.predict(X = xxx.reshape(-1, 1) ), color='red')
        if trendline:
            lbl = "{} = {:.2f} x {} + {:.2f} | R2 = {:.2f}".format(y, p_opt1[0], x , p_opt1[1], r2_1)
            ax[i].plot(xxx, linear(xxx, *p_opt1), 'k--', label=None)
            if trenderror:
                ax[i].plot(xxx, linear(xxx, *p_opt1) + 1.96* np.array([std_linear(xx, p_opt1, p_cov1) for xx in xxx]), 'k--', alpha=0.7)
                ax[i].plot(xxx, linear(xxx, *p_opt1) - 1.96* np.array([std_linear(xx, p_opt1, p_cov1) for xx in xxx]), 'k--', alpha=0.7)
            print("Linear regression: "+lbl)
            if wSigma:
                lbl = "WithSigma: {} = {:.2f} x {} + {:.2f} | R2 = {:.2f}".format(y, p_opt2[0], x , p_opt2[1], r2_2)
                ax[i].plot(xxx, linear(xxx, *p_opt2), 'b-', label=lbl)
                if trenderror:
                    ax[i].plot(xxx, linear(xxx, *p_opt2) + 1.96* np.array([std_linear(xx, p_opt2, p_cov2) for xx in xxx]), 'b--', alpha=0.5)
                    ax[i].plot(xxx, linear(xxx, *p_opt2) - 1.96* np.array([std_linear(xx, p_opt2, p_cov2) for xx in xxx]), 'b--', alpha=0.5)
                
        ax[i].set_ylim([0-5,100+5])
        ax[i].set_xlabel(x, fontsize=14)
        ax[i].set_ylabel(y, fontsize=14)
        
        #ax[i].annotate(text="{} = {:.2f} x {} + {:.2f}\nR2 = {:.2f}".format(y, sl, x ,it, r2_regHC), xy= (0.5, 0.9), xycoords='axes fraction', fontsize=10)
        
        plt.legend()
        fig.suptitle(name, fontsize=15)
        fig.tight_layout()

        #tsp = str(datetime.now().strftime("%Y-%m-%d_%H-%M"))    
        #plt.savefig('data_analysed/img-fid/'+name+'_'+tsp+'.png',dpi=300)
        #plt.show()
        return fig, ax 
    
    if interactive:
        #plt.clf()
        # we assume an ax was provided
        i=0
        #[l.remove() for l in ax[i].lines]
        ax[i].cla()

        ax[i].errorbar(vXa.flatten(), vYa.flatten(), yerr=sigvYa, fmt="none",
                          barsabove=True, elinewidth=0.5, ecolor='k', solid_capstyle='butt', capthick=0.5, capsize=2)
        ax[i].plot(np.array([np.min(vXa),np.max(vXa)]), reg_hc.predict(X = np.array([np.min(vXa),np.max(vXa)]).reshape(-1, 1) ), color='red')
        
        ax[i].plot(xxx, linear(np.array([np.min(vXa),np.max(vXa)]), *p_opt1), 'k--')
        ax[i].plot(xxx, linear(xxx, *p_opt1) + 1.96* np.array([std_linear(xx, p_opt1, p_cov1) for xx in xxx]), 'k--', alpha=0.7)
        ax[i].plot(xxx, linear(xxx, *p_opt1) - 1.96* np.array([std_linear(xx, p_opt1, p_cov1) for xx in xxx]), 'k--', alpha=0.7)
        
        #sc = ax[i].scatter(vXa, vYa, s=20, marker="o")
        sc = sns.scatterplot(data=df, x=x, y=y,
                        style='PyrolysisClass', hue='BiomassClass', palette=sns.color_palette("Set2", as_cmap=False)[:len(df['BiomassClass'].unique())], 
                        ax = ax[i], s=80)
        
        ax[i].set_ylim([0-5,100+5])
        ax[i].set_xlabel(x, fontsize=14)
        ax[i].set_ylabel(y, fontsize=14)
        ax[i].annotate(text="{} = {:.2f} x {} + {:.2f}\nR2 = {:.2f}".format(y, sl, x ,it, r2_regHC), 
                       xy= (0.5, 0.9), xycoords='axes fraction', fontsize=10)
        
        # https://towardsdatascience.com/tooltips-with-pythons-matplotlib-dcd8db758846
        annot = ax[i].annotate("...", xy=(0,0), xytext=(20,20),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(True)

        def update_annot(ind):

            pos = sc.get_offsets()[ind["ind"][0]]
            annot.xy = pos
            text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
            annot.set_text(text)
            annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
            annot.get_bbox_patch().set_alpha(0.4)


        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax[i]:
                cont, ind = sc.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    ifig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        ifig.canvas.draw_idle()

        ifig.canvas.mpl_connect("button_press_event", hover)


## Linear Correlations - timeseries  
def _align_dataframes(fitdata, metadata, outliers):
    '''
    fitdata, sub_metadata, obs = _align_dataframes(fitdata, metadata, outliers)
    '''
    fitdata = fitdata.copy(deep=True)
    obs = list(fitdata['ID_obs'].unique())
    # exclude identified outliers specified in *outliers*
    if outliers is not None:
        obs = set(obs)-set(outliers)
        
    fitdata = fitdata.set_index('ID_obs') # change index to ID_obs, so that series have the same index
    
    fitdata = fitdata[fitdata.index.isin(obs)].copy(deep=True) # filter metadata for the observations in the model minus outliers
    sub_metadata = metadata[metadata.index.isin(obs)] # filter metadata for the observations in the model
    sub_metadata = sub_metadata.replace('na',np.NaN) # replaces 'na' strings by NaN, for correct handling 
    
    # Make sure fitdata and sub_metadata have same index elements & order
    fitdata = fitdata.sort_index() # sort index, to be sure it matches with metadata
    sub_metadata = sub_metadata.sort_index() # sort index, to be sure it matches with metadata
    if not fitdata.index.equals(sub_metadata.index):
        # True if “other” is an Index and it has the same elements and order as the calling index; False otherwise.
        raise Exception("Warning: the indexes of fitdata & sub_metadata do not have same elments or same order")

    return fitdata, sub_metadata, obs

def _retrieveX(x, sub_metadata):
    '''
    vX, vXnan = _retrieveX(x, sub_metadata)
    '''
    acc_x = ['H/C_all'] # special values of x that are accepted (e.g. H/C_all > using H/C_org if available, but bridging with H/C_tot if not available)
    
    if len(x) == 1:
        x = x[0]
        if x not in sub_metadata.columns and x not in acc_x: # verify x is in metadata.columns
            raise Exception("x = `{}` must be a column in `metadata` or an accepted value (check documentation)".format(x))
        elif x in sub_metadata.columns:
            vX = sub_metadata[x]
        else: # x is a special value, swith and handle separately
            if x == 'H/C_all':
                vX = sub_metadata['H/C_org'].fillna(sub_metadata['H/C_tot'], inplace=False) 
        # Counting NaNs in vX
        vXnan = list(vX[vX.isna()].index)
        
    elif len(x) > 1:
       # if all(item in sub_metadata.columns for item in x):
       #     # yes, all items of x are in sub_metadata.columns
       #     vX = sub_metadata[x]
       # elif any(item in acc_x for item in x):
        vX = pd.DataFrame({})
        for item in x:
            if item in acc_x:
                if item == 'H/C_all':
                    vvX = sub_metadata['H/C_org'].fillna(sub_metadata['H/C_tot'], inplace=False) 
                    vX[item] = vvX
            elif item in sub_metadata.columns:
                    vvX = sub_metadata[item]
                    vX[item] = vvX
            else:
                raise Warning("x item `{}` not selected, because it is not in `metadata` or an accepted value (check documentation)".format(item))

        #counting NaNs in vX
        vXnan = list(vX[vX.isna().any(axis=1)].index)
    
    return vX, vXnan  

def _calculateY(fitdata, sub_metadata, Ts, Th):
    '''
    vY, sigvY, vYnan = _calculateY(fitdata, sub_metadata, Ts, Th)
    '''
    # some constants 
    map_model_function = {"singleExp":singleExp,"singleExp_u":singleExp_u,"doubleExp":doubleExp,"doubleExp_u":doubleExp_u,"tripleExp":tripleExp,"tripleExp_u":tripleExp_u,"powerModel":powerModel,}            
    
    # initialise lists to save calcs
    vY = []
    sigvY =[]
    for obs, row in fitdata.iterrows():
        Ts_ini = sub_metadata.loc[obs, 'IncubationTemperature']
        FIT_func = row['model']
        params = row[['k', 'c', 'k1', 'k2', 'c1', 'c2', 'k3', 'c0', 'b', 'm', 'c3']]
        # std is in absolute value (other named .1 => relative value, in %)
        params_std = row[['std_k', 'std_c','std_k1', 'std_k2', 'std_c1', 'std_c2', 'std_k3', 'std_c3', 'std_c0', 'std_b', 'std_m']]
        params_corr = ['corr(k, c)', 'corr(k1, k2)', 'corr(k1, c1)', 'corr(k2, c1)', 'corr(k1, c2)', 'corr(k2, c2)', 'corr(c1, c2)', 'corr(k1, k3)',
                        'corr(k2, k3)', 'corr(k3, c1)', 'corr(k3, c2)', 'corr(k1, c3)', 'corr(k2, c3)', 'corr(k3, c3)', 'corr(c1, c3)', 'corr(c2, c3)','corr(c0, b)', 'corr(c0, m)', 'corr(b, m)']
        params_cov = row[['cov(k, c)', 'cov(k1, k2)', 'cov(k1, c1)', 'cov(k2, c1)', 'cov(k1, c2)', 'cov(k2, c2)',
                       'cov(c1, c2)', 'cov(k1, k3)', 'cov(k2, k3)', 'cov(k3, c1)','cov(k3, c2)', 'cov(k1, c3)', 'cov(k2, c3)', 'cov(k3, c3)','cov(c1, c3)', 'cov(c2, c3)', 'cov(c0, b)', 'cov(c0, m)', 'cov(b, m)']]
        p_opt_names, p_opt, p_cov = rebuild_cov2(map_model_function[FIT_func], params, params_std, params_cov)
        BC, sBC = apply_Q10_bis(Th, Ts, Ts_ini, FIT_func, params, p_opt, p_cov)
        vY.append(BC)
        sigvY.append(sBC[0]) # [0] because functions returns a time series

    vY = pd.Series(data=vY, index=fitdata.index)
    sigvY = pd.Series(data=sigvY, index=fitdata.index)
    vYnan = list(vY[vY.isna()].index)
    return vY, sigvY, vYnan

def linear(x, a, b):
    return a*x+b

def std_linear(x, p_opt, p_cov):
    a = p_opt[0]
    b = p_opt[1]
    gX = np.array([x, 1])
    std_f2 = gX.T @ p_cov @ gX
    return np.sqrt(std_f2)

def calc_correlation_timeseries(fitdata, metadata,
                                x = ['H/C_all'],
                                outliers=None,
                                Ts=[10, 15, 20],
                                Th=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200],
                               ):
    '''
    Function can be used to calculate timesesries of BC_TH_TS, based on a linear fit between BC_TH_TS and H/C_all, for any TH and TS. 

    Inputs:
    - Ts: a list of values, or a single value provided as a list- soil temperatures to consider
    - Th: a list of values, or a single value provided as a list- time points to calculate
    - x : a list of values, or a single value provided as a list - variable in metadata to use for correlation (single variable or multiple variable)
    - outliers: a list of ID_obs to exclude from analysis

    Outputs:
    - vX: vector of X variable, for selected observations where data available
    - vXnan: list of ID_obs where X variable is NaN
    - megaY: dictionnary of calculated Y (BC_TH_TS), where keys are (TS, TH) tuples, and values are the corresponding BC_TH_TS for each observations
    - megaYnan: list of ID_obs where Y variable is NaN
    - megaCoeffs: for each (TS, TH), coefficients from the linear regression (slope, intercept)
    - megaSCoeffs: for each (TS, TH), covariance matrix from the linear regression

    Usage: \n
        vX, vXnan, megaY, megaYnan, megaCoeffs, megaSCoeffs = calc_correlation_timeseries(fitdata, metadata,
                                x = ['H/C_all'], outliers=None, Ts=[10, 15, 20],
                                Th=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200],
                               )
    '''
    
    # data preparation
    fitdata, sub_metadata, obs = _align_dataframes(fitdata, metadata, outliers)
    # get X columns
    vX, vXnan = _retrieveX(x, sub_metadata)
    
    # calculate Y columns: i.e. calculate BC_Th_Ts & sigmaBC_Th_Ts, and save for range of Ts and Th 
    megaY = {} 
    megaSY = {}
    megaYnan = {}
    for ts in Ts:
        for th in Th:
            megaY[(ts, th)], megaSY[(ts, th)], megaYnan[(ts, th)] = _calculateY(fitdata, sub_metadata, ts, th)
            
    # perform correlation, for each Ts and Th; and save coefficients & their sigma
    megaCoeffs = {}
    megaSCoeffs = {}        
    for ts in Ts:
        for th in Th:
            # Remove any mutual observation with nan values
            vXYnan = set(vXnan + megaYnan[(ts, th)])
            cvX = vX.drop(vXYnan, axis=0).sort_index(axis=0)
            cvY = megaY[(ts, th)].drop(vXYnan, axis=0).sort_index(axis=0)
            sigvY = megaSY[(ts, th)].drop(vXYnan, axis=0).sort_index(axis=0) # removing irrelevant observations, as in vX, vY
            sigvY = sigvY.fillna(0) # remainin nan, where uncertainty could not be calculated for some reason, replaced by 0
            # Do the fit
            p_opt1, p_cov1 = curve_fit(f=linear, xdata=cvX, ydata=cvY, p0=None, method='lm', sigma=None, absolute_sigma=False)
            # Save the fit
            megaCoeffs[(ts, th)] = p_opt1
            megaSCoeffs[(ts, th)] = p_cov1
            
    # return array of BC_Th_Ts, sigmaBC_Th_Ts, coefficients & sigma, for each Ts & Th
    # second function: plot all the lines. 
    return vX, vXnan, megaY, megaYnan, megaCoeffs, megaSCoeffs

def _rebuild_timeseries(xx, ts, Th, megaCoeffs):
    BC = []
    for th in Th:
        # get the coefficients
        p_opt = megaCoeffs[(ts, th)]
        BC.append(max(0, min(100, linear(xx, *p_opt))) )
    return BC
    
def plot_correlation_timeseries(fitdata, metadata, x, Ts, Th,
                                vX, vXnan, megaY, megaYnan, megaCoeffs, megaSCoeffs,
                                fig=None, ax=None, TH_lim=1000,
                               ):
    '''
    Plot the results computed via `calc_correlation_timeseries`.
    
    Inputs:
    - x = values of X that need to be plotted, as a list // if multiple variables; should be a dictionnary of key:[values]

    Outputs:
    - fig, ax

    '''
    if (fig, ax) == (None, None):
        fig, ax = plt.subplots(1,1, figsize=(9,5))
    ax.cla()
    # reconstruct the lines, for each H/C_all values, at given Ts
    for xx in x: 
        for ts in Ts:
            BC = _rebuild_timeseries(xx, ts, Th, megaCoeffs)
            ax.plot(Th, BC, label=r"H/C: {}".format(xx))
    #ax.plot(Th, Th)
    
    # styling    
    ax.set_ylim([0-5,100+5])
    ax.set_ylabel('Biochar C remaining (%)', fontsize=14)
    ax.set_xlim([0-10, TH_lim+10])
    ax.set_xlabel('Time (years)', fontsize=14)
    
    ax.legend(fontsize=8, loc=1)
    name="Modelled timeseries as a function of H/C ratio and soil temperature".format(x)
    fig.suptitle(name, fontsize=15)

    return fig, ax

def pick_model(k, df, sets, X, Y):
    '''
    Pick the a random forest model `k` from `df` and re-trains it with its optimal parameters, on the data set (X,Y).

    - k: row index of the model in the dataframe df, corresponding also to set index in dictionary sets
    - df: df containing the model optimal parameters, as well as r2, r2a, features
    - sets: dictionnary of sets 
    - X, Y: the training data/target

    '''
    mp = { k.removeprefix("randomforestregressor__"):v for k,v in df.loc[k,:]['model_params'].items()}

    # feature selection, normalisation, one-hot encoding
    fX = X[list(sets[k])]
    fY = Y

    numerical_columns_selector = selector(dtype_exclude=object)
    categorical_columns_selector = selector(dtype_include=object)
    numerical_columns = numerical_columns_selector(fX)
    categorical_columns = categorical_columns_selector(fX)
    categorical_preprocessor = OneHotEncoder(handle_unknown="ignore") # Here, could do a LabelEncoder() as well
    numerical_preprocessor = RobustScaler() #StandardScaler() # RobustScaler() can also do robustscaler, dealing with outliers differently
    preprocessor = ColumnTransformer([
        ('one-hot-encoder', categorical_preprocessor, categorical_columns),
        ('standard_scaler', numerical_preprocessor, numerical_columns)])

    pipeRF = make_pipeline(preprocessor, RandomForestRegressor(**mp))
    pipeRF.fit(fX, fY)
    return pipeRF