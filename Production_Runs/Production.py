### Optimization of RBM ###
from __future__ import division, print_function, absolute_import
import os
import sys
sys.path.append('/glade/work/dblaskey/RBM_opt_code/src')
import MOASMO
import Gen_Inputs
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

# preprocess RBM
def preprocess_rbm(folder, outlet_id):
    with open('/glade/u/home/dblaskey/RBM/Production_Runs/prepare_rbm_opt_%s_%s.cfg'%(folder,outlet_id), 'w') as f:
        f.write('\n'.join(pre_lines))
        f.write("\n")
        f.write("network_csv: /glade/scratch/dblaskey/RBM/Input/%s_network.csv\n"%outlet_id)
        f.write("fdr_nc: /glade/scratch/dblaskey/RBM/Input/ntopo_MERIT_Hydro_v0.AK_subbasin.%s.nc\n"%outlet_id)
        f.write("flow_output_nc: /glade/scratch/dblaskey/RBM/Output/mizuRoute_Output/AK_Rivers_%s.h.\n"%outlet_id)
        f.write("grid_weight: /glade/scratch/dblaskey/RBM/Input/spatialweights_AK_to_merit-%s.nc\n"%outlet_id)
        f.write("input_id: %s\n"%outlet_id)
        f.write("site_loc: %s\n"%local)
        f.write("optimization_param: /glade/u/home/dblaskey/RBM/Production_Runs/runs.csv\n")
        f.write("rbm_energy_file: /glade/scratch/dblaskey/RBM/RBM_Input/%s_energy\n"%outlet_id)
        
    config_file ='/glade/u/home/dblaskey/RBM/Production_Runs/prepare_rbm_opt_%s_%s.cfg'%(folder,outlet_id)
    Gen_Inputs.generate_RBM_inputs(config_file, folder)

# Execute RBM
def run_rbm(folder, outlet_id):
    file_path = '/glade/scratch/dblaskey/RBM/Output/%s/%s_2021.temp'%(folder, outlet_id)

    if os.path.isfile(file_path) == True:
        print("Skip ", outlet_id, " in ", folder)

    else:
        print("Running ", outlet_id, " in ", folder)
        # Call function to run RBM
        args ='/glade/work/dblaskey/RBM/src/rbm10_mizu /glade/scratch/dblaskey/RBM/RBM_Input/%s/%s /glade/scratch/dblaskey/RBM/Output/%s/%s'%(folder, outlet_id, folder, outlet_id) 
        subprocess.run(args, shell=True)

# Post-process RBM by converting text files to NetCDF     
def postprocess_rbm(folder, outlet_id):
    # Specify path
    file_path = '/glade/scratch/dblaskey/RBM/Output/%s/%s_2021.nc'%(folder, outlet_id)
    # Convert output to NetCDF file                  
    if os.path.isfile(file_path) == True:
        print("Skip making nc for ", outlet_id, " in ", folder)
    else:
        print("Running text2nc for ", outlet_id, " in ", folder)

        for year in years:
            with open('/glade/u/home/dblaskey/RBM/Production_Runs/postprocess_rbm_opt_%s_%s.cfg'%(folder,outlet_id), 'w') as f:
                f.write('[INPUT]\n')
                f.write('spat_file: /glade/scratch/dblaskey/RBM/Output/%s/%s.Spat\n'%(folder, outlet_id))
                f.write("temp_file: /glade/scratch/dblaskey/RBM/Output/%s/%s_%s.temp\n"%(folder, outlet_id, year))
                f.write("[PARAM]\n")
                f.write("sim_year: %s\n"%year)
                f.write("[OUTPUT]\n")
                f.write("rbm_nc_file: /glade/scratch/dblaskey/RBM/Output/%s/%s_%s.nc"%(folder, outlet_id, year))
            
            config_file ='/glade/u/home/dblaskey/RBM/Production_Runs/postprocess_rbm_opt_%s_%s.cfg'%(folder, outlet_id)
            convert_rbm_to_nc.rbm_to_nc(config_file)

#### Begin Code #####
df_sites = pd.read_csv('/glade/u/home/dblaskey/RBM/Production_Runs/AK_outlets.csv', index_col=0)
outlets =  [81030789, 81035794, 81036242, 81036460, 81036473] #np.unique(df_sites.COMID.values)

start_date = '1990-01-01'
end_date = '2021-09-30'
years = range(1990,2022)
time_series = pd.date_range(start_date, end_date)

pre_lines = ["[INPUT]", 
"output_energy_nc: /glade/scratch/tcraig/archive/NNA.4km.hERA5.1989.003.pp/lnd/hist/NNA.4km.hERA5.1989.003.clm2.hrbm2.", "","[PARAM]", "start_date: 19900101", "end_date: 20210930", "seg_minflow : 4", "sd_time_step: 1", "", "[RBM_OPTIONS]", "min_flow : 5", "", "[OPT]"]

location="pe_basin"

psets_df = pd.read_csv('/glade/u/home/dblaskey/RBM/Production_Runs/runs.csv')
psets_df.rename(columns={ psets_df.columns[0]: "Name" }, inplace = True)

# create a list of arguments for the functions and build folders
args_list = []
for folder in np.unique(psets_df.Name.values):
    # Build Input Folders
    path = '/glade/scratch/dblaskey/RBM/RBM_Input/%s/'%folder
    if os.path.exists(path) == True:
        print("Skip Creating ", folder)            
    else:             
        print("Creating ", folder) 
        os.mkdir(path) 

    # Build Output Folders
    path = '/glade/scratch/dblaskey/RBM/Output/%s/'%folder
    if os.path.exists(path) == True:
        print("Skip Creating ", folder)            
    else:
        print("Creating ", folder)
        os.mkdir(path)

    for outlet in outlets:
        args_list.append((folder, outlet))

for local in np.unique(df_sites.Location.values):
    
    df_site_local = df_sites[df_sites['Location']==local]
    psets_df_local = psets_df[psets_df["Location"] == local]
    
    ### Build argument list ###
    args_list_prep = []
    for folder_local in psets_df_local.Name.values:
        for outlet_local in [81030789, 81035794, 81036242, 81036460, 81036473]: #np.unique(df_site_local.COMID.values):
            args_list_prep.append((folder_local, outlet_local))
    
    #### Preprocess RBM #### 
    # set up the multiprocessing pool
    pool = multiprocessing.Pool()

    # run the function in parallel using the multiprocessing pool
    results = pool.starmap(preprocess_rbm, args_list_prep)

    # close the multiprocessing pool
    pool.close() 

    # Clean up    
    sys_args ='rm -r /glade/u/home/dblaskey/RBM/Production_Runs/prepare_rbm_opt_*.cfg'
    subprocess.run(sys_args, shell=True)

    print("Done with pre-procession")

#### Run RBM #### 
# set up the multiprocessing pool
pool = multiprocessing.Pool()

# run the function in parallel using the multiprocessing pool
results = pool.starmap(run_rbm, args_list)

# close the multiprocessing pool
pool.close()   

print("Done with RBM runs")

# ## Post process RBM ####

# set up the multiprocessing pool
pool = multiprocessing.Pool()

# run the function in parallel using the multiprocessing pool
results = pool.starmap(postprocess_rbm, args_list)

# close the multiprocessing pool
pool.close()   

# Clean Up
args ='rm -r /glade/u/home/dblaskey/RBM/Production_Runs/postprocess_rbm_opt_*.cfg'
subprocess.run(args, shell=True)

print("Done with text to netCDF conversion")
