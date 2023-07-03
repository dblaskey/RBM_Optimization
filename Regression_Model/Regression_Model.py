### Optimization of RBM ###
from __future__ import division, print_function, absolute_import
import os
import sys
sys.path.append('/glade/work/dblaskey/RBM_opt_code/src')
import MOASMO
import Gen_Val_Inputs
import convert_rbm_to_nc
from convert_rbm_to_nc import read_config
import numpy as np
import matplotlib.pyplot as plt
import util
import _pickle as cPickle
import pandas as pd
import sampling
import xarray as xr
import seaborn as sn
import GLP
import gp
import NSGA2
import pickle
import smt
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
import random
import re
from sklearn.metrics import mean_squared_error
import cmocean
import hydroeval as he
import glob
import subprocess
import multiprocessing
import math
import warnings
from scipy.ndimage import gaussian_filter1d

#### Begin Code #####
df_sites = pd.read_csv('/glade/u/home/dblaskey/RBM/Regression_Model/Opt_Basins.csv', index_col=0)
df_sites = df_sites[df_sites['type']=="Opt"]
outlets = np.unique(df_sites.outlet_comid.values)

# Read in Observed Data
temp_data = pd.read_csv('/glade/scratch/dblaskey/RBM/temperature_gages.csv', index_col=0)
df = pd.merge(temp_data, df_sites, on="site_no")

psets_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs_0.csv')
psets_df.rename(columns={ psets_df.columns[0]: "Name" }, inplace = True)

years = range(2013,2018)
start_date = '2013-01-01'
end_date = '2017-12-31'
time_series = pd.date_range(start_date, end_date)

for folder in np.unique(psets_df.Name.values):
    # Build Output Folders
    path = '/glade/scratch/dblaskey/RBM/Regression/Opt/test/%s/'%folder
    if os.path.exists(path) == False:           
        print("Creating ", folder)
        os.mkdir(path)
        
# Calculate air temperature
def mosheni(air_temp, alpha, beta, mu, gamma):
    return mu + (alpha - mu)/(1 + math.exp(gamma*(beta-air_temp)))

#### Begin Code #####
df_sites = pd.read_csv('/glade/u/home/dblaskey/RBM/Regression_Model/Opt_Basins.csv', index_col=0)
df_sites = df_sites[df_sites['type']=="Opt"]
outlets = np.unique(df_sites.outlet_comid.values)

# Read in Observed Data
temp_data = pd.read_csv('/glade/scratch/dblaskey/RBM/temperature_gages.csv', index_col=0)
df = pd.merge(temp_data, df_sites, on="site_no")

psets_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs_0.csv')
psets_df.rename(columns={ psets_df.columns[0]: "Name" }, inplace = True)

years = range(2013,2018)
start_date = '2013-01-01'
end_date = '2017-12-31'
time_series = pd.date_range(start_date, end_date)

for folder in np.unique(psets_df.Name.values):
    # Build Output Folders
    path = '/glade/scratch/dblaskey/RBM/Regression/Opt/test/%s/'%folder
    if os.path.exists(path) == False:           
        print("Creating ", folder)
        os.mkdir(path)

def cal_regression(folder):
    for outlet in outlets:
        print(folder, outlet)
        
        var_list = psets_df[psets_df['Name'] == folder]

        alpha = var_list.alpha.values[0]
        beta = var_list.beta.values[0]
        mu = var_list.mu.values[0]
        gamma = var_list.gamma.values[0]

        Ordered_reaches_final = pd.read_csv("/glade/scratch/dblaskey/RBM/Input/%s_network.csv" % outlet)
        
        results = []

        for i, comid in enumerate(df_sites[df_sites['outlet_comid'] == outlet]['COMID'].values):
            cell = Ordered_reaches_final[Ordered_reaches_final['hru'] == comid].node.values[0]
            site_no = df_sites[df_sites['outlet_comid'] == outlet].index.values[i]
            
            path_results = '/glade/scratch/dblaskey/RBM/Regression/Opt/test/%s/%s.csv' % (folder, comid)
            if os.path.exists(path_results) == False:
                Water = []
                for year in years:
                    test = pd.read_csv("/glade/scratch/dblaskey/RBM/RBM_Input/%s_energy_%s" % (outlet, year), sep=" ", header=None, names=["cell", "Tair", "vp", "SW", "LW", "Density", "P", "Wind"])
                    test_id = test[test["cell"] == cell]
                    test_id.drop(test_id.tail(1).index, inplace=True)
                    test_id['T_smooth'] = 0.1 * test_id['Tair'] + 0.9 * test_id['Tair'].shift()
                    test_id['T_smooth'] = test_id['T_smooth'].fillna(test_id['Tair'])

                    # Calculate air temperature
                    def mosheni(air_temp, alpha, beta, mu, gamma):
                        return mu + (alpha - mu)/(1 + math.exp(gamma*(beta-air_temp)))


                    Water.append(test_id.apply(lambda row : mosheni(row['T_smooth'], alpha, beta, mu, gamma), axis=1))

                Water = np.concatenate(Water)                        
                sim_df = pd.DataFrame(Water, index=time_series, columns=['sim'])
                sim_df.index.name = 'Date'

                # Filter to just one gage for observations
                temp_df = df[df['site_no'] == site_no][['X_00010_00003', 'Date']]
                temp_df = temp_df.rename(columns={"X_00010_00003": "obs"})
                temp_df = temp_df.set_index('Date')
                temp_df.index = pd.to_datetime(temp_df.index)

                # Combine simulation and observation datasets
                df_concat = pd.concat([sim_df, temp_df], axis=1)
                df_concat = df_concat.dropna()
                df_concat = df_concat.loc['2013-10-01':'2017-09-30']
                df_concat = df_concat[df_concat.index.month.isin([5, 6, 7, 8, 9])]
                df_concat = df_concat.reset_index()

                df_concat.to_csv('/glade/scratch/dblaskey/RBM/Regression/Opt/test/%s/%s.csv' % (folder, comid))

                result = {
                    'RMSE': he.evaluator(he.rmse, df_concat.sim.values, df_concat.obs.values)[0],
                    'COMID': np.int(comid),
                    'Name': folder,
                    'Outlet': outlet
                }

                results.append(result)
            
            else:
                df_concat = pd.read_csv(path_results)

                result = {
                    'RMSE': he.evaluator(he.rmse, df_concat.sim.values, df_concat.obs.values)[0],
                    'COMID': np.int(comid),
                    'Name': folder,
                    'Outlet': outlet
                }

                results.append(result)
    return results

with warnings.catch_warnings(record=True):
    #### Preprocess RBM #### 
    # set up the multiprocessing pool
    pool = multiprocessing.Pool()

    # run the function in parallel using the multiprocessing pool
    results = pool.map(cal_regression, np.unique(psets_df.Name.values))

    # close the multiprocessing pool
    pool.close() 

    warnings.warn("should not appear")

# Flatten the list of lists into a single list
results = [item for sublist in results for item in sublist]

# Convert the list of dictionaries to a DataFrame
results_df = pd.DataFrame(results)

df2 = pd.merge(results_df, df_sites, on=['COMID'])
df2.to_csv('/glade/u/home/dblaskey/RBM/Regression_Model/LHS_results.csv', index=False)

df_LHS = results_df.groupby('Name').mean("temp_rmse").reset_index().drop(columns=['COMID', 'Outlet'])
df_LHS = df_LHS.merge(psets_df)

df_LHS = df_LHS.drop(["min_velocity", "min_d", "min_w"], axis = 1)

df_LHS.to_csv('/glade/u/home/dblaskey/RBM/Regression_Model/LHS_results_full.csv', index=False)

