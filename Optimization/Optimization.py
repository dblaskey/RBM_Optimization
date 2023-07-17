### Optimization of RBM ###
from __future__ import division, print_function, absolute_import
import os
import sys
sys.path.append('/glade/work/dblaskey/RBM_opt_code/src')
import MOASMO
import Gen_RBM_Inputs
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
from functools import partial

# define the function that will be run in parallel

# preprocess RBM
def preprocess_rbm(folder, outlet_id):
    with open('/glade/u/home/dblaskey/RBM/Optimization/prepare_rbm_opt_%s_%s.cfg'%(folder,outlet_id), 'w') as f:
        f.write('\n'.join(pre_lines))
        f.write("\n")
        f.write("network_csv: /glade/scratch/dblaskey/RBM/Input/%s_network.csv\n"%outlet_id)
        f.write("fdr_nc: /glade/scratch/dblaskey/RBM/Input/ntopo_MERIT_Hydro_v0.AK_subbasin.%s.nc\n"%outlet_id)
        f.write("flow_output_nc: /glade/scratch/dblaskey/RBM/Output/mizuRoute_Output/AK_Rivers_%s.h.\n"%outlet_id)
        f.write("grid_weight: /glade/scratch/dblaskey/RBM/Input/spatialweights_AK_to_merit-%s.nc\n"%outlet_id)
        f.write("input_id: %s\n"%outlet_id)
        f.write("optimization_param: /glade/u/home/dblaskey/RBM/Optimization/Opt_runs_%s.csv\n"%iteration)
        f.write("rbm_energy_file: /glade/scratch/dblaskey/RBM/RBM_Input/Calibration/%s_energy\n"%outlet_id)
        
    config_file ='/glade/u/home/dblaskey/RBM/Optimization/prepare_rbm_opt_%s_%s.cfg'%(folder,outlet_id)
    Gen_RBM_Inputs.generate_RBM_inputs(config_file, folder)

# Execute RBM
def run_rbm(folder, outlet_id):
    file_path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_2017.temp'%(folder, outlet_id)

    if os.path.isfile(file_path) == True:
        print("Skip ", outlet_id, " in ", folder)

    else:
        print("Running ", outlet_id, " in ", folder)
        # Call function to run RBM
        args ='/glade/work/dblaskey/RBM/src/rbm10_mizu /glade/scratch/dblaskey/RBM/RBM_Input/Optimize/%s/%s /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s'%(folder, outlet_id, folder, outlet_id) 
        subprocess.run(args, shell=True)
    
# Post-process RBM by converting text files to NetCDF     
def postprocess_rbm(folder, outlet_id):
    # Specify path
    file_path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_2017.nc'%(folder, outlet_id)
    # Convert output to NetCDF file                  
    if os.path.isfile(file_path) == True:
        print("Skip making nc for ", outlet_id, " in ", folder)
    else:
        print("Running text2nc for ", outlet_id, " in ", folder)

        for year in years:
            with open('/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm_opt_%s_%s.cfg'%(folder,outlet_id), 'w') as f:
                f.write('[INPUT]\n')
                f.write('spat_file: /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s.Spat\n'%(folder, outlet_id))
                f.write("temp_file: /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_%s.temp\n"%(folder, outlet_id, year))
                f.write("[PARAM]\n")
                f.write("sim_year: %s\n"%year)
                f.write("[OUTPUT]\n")
                f.write("rbm_nc_file: /glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_%s.nc"%(folder, outlet_id, year))
            
            config_file ='/glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm_opt_%s_%s.cfg'%(folder, outlet_id)
            convert_rbm_to_nc.rbm_to_nc(config_file)
            
# Calculate optimzation results
def process_folder_outlet(folder, outlet):
    print(folder, outlet)
    gages = df_sites.loc[df_sites['outlet_comid'] == outlet]
    comid = gages.COMID.values
    site_no = gages.index.values

    # Specify path
    ds_rbm = xr.open_mfdataset('/glade/scratch/dblaskey/RBM/Output/Optimize/%s/%s_*.nc'%(folder, outlet))
    reachID_list = np.rint(ds_rbm.hru.values).astype(int)

    # Filter to just one gage for simulations
    seg_sel = np.where(reachID_list == comid)[0][0]
    sim_flow = ds_rbm.T_stream.sel(hru=reachID_list[seg_sel], no_seg=2)
    sim_flow.load()
    sim_df = pd.DataFrame(sim_flow, index=time_series, columns=['sim'])
    sim_df.index.name = 'Date'

    # Filter to just one gage for observations
    temp_df = df[df['site_no'].isin(site_no)][['X_00010_00003', 'Date']]
    temp_df = temp_df.rename(columns={"X_00010_00003": "obs"})
    temp_df = temp_df.set_index('Date')
    temp_df.index = pd.to_datetime(temp_df.index)

    # Combine simulation and observation datasets
    df_concat = pd.concat([sim_df,temp_df],axis=1)
    df_concat = df_concat.dropna()
    df_concat = df_concat.loc['2013-10-01':'2017-09-30']
    df_concat = df_concat[df_concat.index.month.isin([4,5,6,7,8,9,10])]
    df_concat = df_concat.reset_index()

    # return the results
    return [he.evaluator(he.rmse, df_concat.sim.values, df_concat.obs.values)[0],
            he.evaluator(he.pbias, df_concat.sim.values, df_concat.obs.values)[0],
            np.int(comid), folder, outlet]

#### Begin Code #####
param_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Param_range.csv')

df_sites = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_Basins.csv', index_col=0)
df_sites = df_sites[df_sites['type'] == "Opt"]

# Read in Observed Data
temp_data = pd.read_csv('/glade/scratch/dblaskey/RBM/temperature_gages.csv', index_col=0)
df = pd.merge(temp_data, df_sites, on="site_no")

outlets = df_sites.outlet_comid.values

start_date = '2013-01-01'
end_date = '2017-12-31'
years = range(2013,2018)
time_series = pd.date_range(start_date, end_date)

pre_lines = ["[INPUT]", 
"output_energy_nc: /glade/scratch/tcraig/archive/NNA.4km.hERA5.1989.003.pp/lnd/hist/NNA.4km.hERA5.1989.003.clm2.hrbm2.", "","[PARAM]", "start_date: 20130101", "end_date: 20171231", "seg_minflow : 2", "sd_time_step: 1", "", "[RBM_OPTIONS]", "min_flow : 5", "", "[OPT]"]

location="pe_basin"

# Initiate Latin Hyper Cube
if os.path.exists('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs_0.csv') == False:
    nInput = len(param_df)
    pct = 0.2
    pop = 100
    Xinit = None
    Yinit = None
    N_resample = int(pop*pct)
    Ninit = 200 #
    Xinit = sampling.glp(Ninit, nInput)

    xub = param_df['Max_Value'].values
    xlb = param_df['Min_Value'].values

    perturbed_param = []
    for i in range(Ninit):
        perturbed_param.append(Xinit[i,:] * (xub - xlb) + xlb)

    perturbed_param = np.array(perturbed_param)
    iteration=0

    test_id_list=[]
    for id_ in range(Ninit):
        test_id = '%s_%s_%s'%(location, "%i"%(iteration), "%04i"%(id_))
        test_id_list.append(test_id)

    psets_df_0 = pd.DataFrame(perturbed_param, columns=param_df['Var_name'].values, index=test_id_list)
    psets_df_0.to_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs_0.csv')

for iteration in range(12,16):
    if iteration == 0:
        print('Starting on Latin Hyper Cube')
    else:
        print("Starting iteration: ", iteration)
        df_LHS_0 = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/LHS_results_0.csv')

        d = preprocessing.normalize(df_LHS_0.drop(["Name", "temp_rmse"], axis = 1),axis=0,return_norm=True)
        normalization_scalar = d[1]
    
        df_LHS_pre = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/LHS_results_%s.csv'%(iteration - 1))
    
    psets_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs_%s.csv'%iteration)
    psets_df.rename(columns={ psets_df.columns[0]: "Name" }, inplace = True)
    
    # create a list of arguments for the functions and build folders
    args_list = []
    for folder in psets_df.Name.values:
        # Build Input Folders
        path = '/glade/scratch/dblaskey/RBM/RBM_Input/Optimize/%s/'%folder
        if os.path.exists(path) == True:
            print("Skip Creating ", folder)            
        else:             
            print("Creating ", folder) 
            os.mkdir(path) 
        
        # Build Output Folders
        path = '/glade/scratch/dblaskey/RBM/Output/Optimize/%s/'%folder
        if os.path.exists(path) == True:
            print("Skip Creating ", folder)            
        else:
            print("Creating ", folder)
            os.mkdir(path)
        
        for outlet in outlets:
            args_list.append((folder, outlet))
    
    #### Preprocess RBM #### 
    # set up the multiprocessing pool
    pool = multiprocessing.Pool()

    # run the function in parallel using the multiprocessing pool
    results = pool.starmap(preprocess_rbm, args_list)

    # close the multiprocessing pool
    pool.close() 
    
    # Clean up    
    sys_args ='rm -r /glade/u/home/dblaskey/RBM/Optimization/prepare_rbm_opt_*.cfg'
    subprocess.run(sys_args, shell=True)
    
    print("Done with pre-procession iteration: ", iteration)
    
    #### Run RBM #### 
    # set up the multiprocessing pool
    pool = multiprocessing.Pool()

    # run the function in parallel using the multiprocessing pool
    results = pool.starmap(run_rbm, args_list)

    # close the multiprocessing pool
    pool.close()   
    
    print("Done with RBM runs for iteration: ", iteration)
    
    ### Post process RBM ####
    
    # set up the multiprocessing pool
    pool = multiprocessing.Pool()

    # run the function in parallel using the multiprocessing pool
    results = pool.starmap(postprocess_rbm, args_list)

    # close the multiprocessing pool
    pool.close()   
    
    # Clean Up
    args ='rm -r /glade/u/home/dblaskey/RBM/Optimization/postprocess_rbm_opt_*.cfg'
    subprocess.run(args, shell=True)
    
    print("Done with text to netCDF conversion for iteration: ", iteration)

    #### Calculate Optimization Metric ####     
    # set up the multiprocessing pool
    pool = multiprocessing.Pool()

    # run the function in parallel using the multiprocessing pool
    results = pool.starmap(process_folder_outlet, args_list)

    # close the multiprocessing pool
    pool.close()

    # create a dataframe from the results
    df_temp = pd.DataFrame(results, columns=['temp_rmse', 'temp_pbias', 'COMID', 'Name', 'Outlet'])           
    df2 = pd.merge(df_temp, psets_df, on=['Name'])
    df2 = pd.merge(df2, df_sites, on=['COMID'])
    df2.to_csv('/glade/u/home/dblaskey/RBM/Optimization/LHS_results_raw_loc_%s.csv'%iteration, index=False)
    
    df_temp = df_temp.groupby('Name').mean("temp_rmse").reset_index().drop(columns=['COMID', 'temp_pbias', 'Outlet'])
    df_LHS = df_temp.merge(psets_df)
    
    if iteration == 0:
        d = preprocessing.normalize(df_LHS.drop(["Name", "temp_rmse"], axis = 1),axis=0,return_norm=True)
        normalization_scalar = d[1]
    else:
        df_LHS = pd.concat([df_LHS_pre, df_LHS]).reset_index(drop=True)
    
    df_LHS.to_csv('/glade/u/home/dblaskey/RBM/Optimization/LHS_results_%s.csv'%iteration, index=False)
    
    # normalize values
    scaled = df_LHS.drop(["Name", "temp_rmse"], axis = 1)/normalization_scalar
    
    # start training the surrogate models
    nInput, nOutput = len(param_df), 1
    alpha = 1e-4
    lb = 1e-4
    ub = 1e3

    # perform optimization using the surrogate model
    gen = 100
    crossover_rate = 0.9
    mu = 20
    mum = 20
    N_resample = 20
    leng_lb = 1e-4
    leng_ub = 1e3
    nu = 1.5
    pop = 100
    
    # start training the surrogate models
    x = scaled.values
    y = df_LHS["temp_rmse"].values

    xlb_single_value_scaled = param_df['Min_Value']/normalization_scalar
    xub_single_value_scaled = param_df['Max_Value']/normalization_scalar

    sm = gp.GPR_Matern(x, y, nInput, nOutput, x.shape[0], xlb_single_value_scaled, xub_single_value_scaled, alpha=alpha, leng_sb=[leng_lb,leng_ub], nu=nu)

    bestx_sm, besty_sm, x_sm, y_sm = \
        NSGA2.optimization(sm, nInput, nOutput, xlb_single_value_scaled.values, xub_single_value_scaled.values, \
                           pop, gen, crossover_rate, mu, mum)
    D = NSGA2.crowding_distance(besty_sm)
    idxr = D.argsort()[::-1][:N_resample]
    x_resample = bestx_sm[idxr,:]
    y_resample = np.zeros((N_resample,nOutput))

    # create test id
    test_id_list=[]
    for id_ in range(N_resample):
        test_id = '%s_%s_%s'%(location, "%i"%(iteration+1), "%04i"%(id_))
        test_id_list.append(test_id)
    psets_df = pd.DataFrame(x_resample, columns=param_df['Var_name'].values, index=test_id_list)

    psets_df = psets_df*normalization_scalar

    psets_df.to_csv('/glade/u/home/dblaskey/RBM/Optimization/Opt_runs_%s.csv'%(iteration+1))
    print("Done with iteration: ", iteration)