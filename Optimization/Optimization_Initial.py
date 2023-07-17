### Optimization of RBM ###
from __future__ import division, print_function, absolute_import
import os
import sys
sys.path.append('/glade/work/dblaskey/RBM_opt_code/src')
import MOASMO
import Gen_Opt_Inputs
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
from functools import partial

# task to execute in another process
def preprocess_rbm(outlet_id):
    with open('/glade/u/home/dblaskey/RBM/Optimization/prepare_rbm_opt_%s.cfg'%outlet_id, 'w') as f:
        f.write('\n'.join(pre_lines))
        f.write("\n")
        f.write("fdr_nc: /glade/scratch/dblaskey/RBM/Input/ntopo_MERIT_Hydro_v0.AK_subbasin.%s.nc\n"%outlet_id)
        f.write("flow_output_nc: /glade/scratch/dblaskey/RBM/Output/mizuRoute_Output/AK_Rivers_%s.h.\n"%outlet_id)
        f.write("grid_weight: /glade/scratch/dblaskey/RBM/Input/spatialweights_AK_to_merit-%s.nc\n"%outlet_id)
        f.write("input_id: %s\n"%outlet_id)
        f.write("optimization_param: /glade/u/home/dblaskey/RBM/Optimization/Opt_runs.csv\n")
        f.write("rbm_energy_file: /glade/scratch/dblaskey/RBM/RBM_Input/Calibration/%s_energy\n"%outlet_id)
        
    config_file ='/glade/u/home/dblaskey/RBM/Optimization/prepare_rbm_opt_%s.cfg'%outlet_id 
    Gen_Opt_Inputs.generate_RBM_inputs(config_file)

# Execute RBM
def run_rbm(outlet_id):
    file_path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_2017.temp'%(folder, outlet_id)

    if os.path.isfile(file_path) == True:
        print("Skip ", outlet_id, " in ", folder)

    else:
        print("Running ", outlet_id, " in ", folder)
        # Call function to run RBM
        args ='/glade/work/dblaskey/RBM/src/rbm10_mizu /glade/scratch/dblaskey/RBM/RBM_Input/Optimize/%s/%s /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s'%(folder, outlet_id, folder, outlet_id) 
        subprocess.run(args, shell=True)

# Post-process RBM by converting text files to NetCDF     
def postprocess_rbm(outlet_id):
    # Specify path
    file_path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_2017.nc'%(folder, outlet_id)
    # Convert output to NetCDF file                  
    if os.path.isfile(file_path) == True:
        print("Skip making nc for ", outlet_id, " in ", folder)
    else:
        print("Running text2nc for ", outlet_id, " in ", folder)

        for year in years:
            with open('/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm_opt_%s.cfg'%outlet_id, 'w') as f:
                f.write('[INPUT]\n')
                f.write('spat_file: /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s.Spat\n'%(folder, outlet_id))
                f.write("temp_file: /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_%s.temp\n"%(folder, outlet_id, year))
                f.write("[PARAM]\n")
                f.write("sim_year: %s\n"%year)
                f.write("[OUTPUT]\n")
                f.write("rbm_nc_file: /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_%s.nc"%(folder, outlet_id, year))
            # Call function
            args ='python /glade/u/home/dblaskey/RBM/Functions/convert_rbm_to_nc.py /glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm_opt_%s.cfg'%outlet_id 
            subprocess.run(args, shell=True)

# Calculate Optimization Metrics
def calculate_metrics(folder):
    rmse_list = []
    pbias_list = []
    print(folder)

    for outlet in outlets:
        # Read in Simulated Data
        gages = df_sites.loc[df_sites['outlet_comid'] == outlet]
        gages = gages[gages['obs_no_year.y'].notna()]
        comid = gages.COMID.values
        site_no = gages.index.values

        # Specify path
        ds_rbm = xr.open_mfdataset('%s/%s_*.nc'%(folder, outlet))
        reachID_list = np.rint(ds_rbm.hru.values).astype(int)

        for i in range(len(gages)):

            # Filter to just one gage for simulations
            seg_sel = np.where(reachID_list == comid[i])[0][0]
            sim_flow = ds_rbm.T_stream.sel(hru=reachID_list[seg_sel], no_seg=2)
            sim_flow.load()
            sim_df = pd.DataFrame(sim_flow, index=time_series, columns=['sim'])
            sim_df.index.name = 'Date'

            # Filter to just one gage for observations
            temp_df = df[df['site_no'] == site_no[i]][['X_00010_00003', 'Date']]
            temp_df = temp_df.rename(columns={"X_00010_00003": "obs"})
            temp_df = temp_df.set_index('Date')
            temp_df.index = pd.to_datetime(temp_df.index)

            # Combine simulation and observation datasets
            df_concat = pd.concat([sim_df,temp_df],axis=1)
            df_concat = df_concat.dropna()
            df_concat = df_concat.loc['2013-10-01':'2017-09-30']
            df_concat = df_concat[df_concat.index.month.isin([4,5,6,7,8,9,10])]
            df_concat = df_concat.reset_index()

            rmse_list.append(he.evaluator(he.rmse, df_concat.sim.values, df_concat.obs.values)[0])
            pbias_list.append(he.evaluator(he.pbias, df_concat.sim.values, df_concat.obs.values)[0])

    df_temp.at[folder,'temp_rmse'] = np.array(rmse_list).mean()
    df_temp.at[folder,'temp_pbias'] = np.array(pbias_list).mean()

    return np.array(rmse_list).mean(), np.array(pbias_list).mean()

param_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Param_range.csv')

df_sites = pd.read_csv('/glade/scratch/dblaskey/RBM/calibration_points.csv', index_col=0)

# Read in Observed Data
temp_data = pd.read_csv('/glade/scratch/dblaskey/RBM/temperature_gages.csv', index_col=0)

df = pd.merge(temp_data, df_sites, on="site_no")
df = df[df['type'] != "Discharge Gage Only"]
df = df.drop(columns=['obs_no_year.x', 'obs_no_year.y'])

outlets = [81015621, 81015538]

start_date = '2013-01-01'
end_date = '2017-12-31'
years = range(2013,2018)
time_series = pd.date_range(start_date, end_date)

location="pe_basin"

pre_lines = ["[INPUT]", 
"output_energy_nc: /glade/scratch/tcraig/archive/NNA.4km.hERA5.1989.003.pp/lnd/hist/NNA.4km.hERA5.1989.003.clm2.hrbm2.", "","[PARAM]", "start_date: 20130101", "end_date: 20171231", "seg_minflow : 2", "sd_time_step: 1", "", "[RBM_OPTIONS]", "min_flow : 5", "min_velocity: 0.75", "", "[OPT]"]   
    
psets_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs.csv')
psets_df.rename(columns={ psets_df.columns[0]: "Name" }, inplace = True)
    
#### Pre-process RBM ####
if __name__ == '__main__':
    with multiprocessing.Pool() as pool:
        pool.map(preprocess_rbm, outlets)
    print("Done with pre-procession")
    
#### Run RBM ####  
for i in range(len(psets_df)):  
    folder = psets_df.iloc[i,0]

    # Specify path
    path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/'%folder

    # Check whether the specified
    # path exists or not
    if os.path.exists(path) == True:
        print("Skip Creating ", folder)            
    else:
        print("Creating ", folder)
        os.mkdir(path)
         
    # Execute RBM  
    if __name__ == '__main__':
        with multiprocessing.Pool() as pool:
            pool.map(run_rbm, outlets)
    
    # Change text files to NetCDF files
    if __name__ == '__main__':
        with multiprocessing.Pool() as pool:
            pool.map(postprocess_rbm, outlets)
        
print("Done with RBM runs")

#### Calculate Optimization Metric ####
    
#folder_list = glob.glob('/glade/scratch/dblaskey/RBM/Output/Optimize/*')

# Create temp dataframe
#df_temp = pd.DataFrame(columns = ['temp_rmse', 'temp_pbias'], index = folder_list)

#if __name__ == '__main__':
#    with multiprocessing.Pool() as pool:
#        results=pool.map(calculate_metrics, folder_list)
#
#   for i, folder in enumerate(folder_list):
#        df_temp.at[folder, 'temp_rmse'] = results[i][0]
#        df_temp.at[folder, 'temp_pbias'] = results[i][1]
#               
#df_temp = df_temp.rename_axis('folder').reset_index()
#df_temp["Name"]= df_temp.folder.str.split("/", expand = True)[7]
#
#df_LHS = df_temp.merge(psets_df)
#    
#df_LHS.to_csv('/glade/u/home/dblaskey/RBM/Optimization/LHS_results_0.csv', index=False)