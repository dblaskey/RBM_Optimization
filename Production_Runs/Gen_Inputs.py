#!/usr/local/anaconda/bin/python

# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#
# This script creates the network, energy and flow files 
# needed to run RBM
#
# Author: Dylan Blaskey
# Created: 5/26/22
# Last modified: 
#
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

import numpy as np
import sys
import pandas as pd
import xarray as xr
import gc
import timeit
import glob
import subprocess
import os

#===============define functions===============#
def read_config(config_file, default_config=None):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """

    from netCDF4 import Dataset
    try:
        from cyordereddict import OrderedDict
    except:
        from collections import OrderedDict
    try:
        from configparser import SafeConfigParser
    except:
        from ConfigParser import SafeConfigParser
    import configobj

    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = OrderedDict()
    for section in sections:
        options = config.options(section)
        dict2 = OrderedDict()
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2

    if default_config is not None:
        for name, section in dict1.items():
            if name in default_config.keys():
                for option, key in default_config[name].items():
                    if option not in section.keys():
                        dict1[name][option] = key

    return dict1
# -------------------------------------------------------------------- #

def config_type(value):
    """
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        else:
            try:
                return int(value)
            except:
                pass
            try:
                return float(value)
            except:
                return value
    else:
        try:
            return list(map(int, val_list))
        except:
            pass
        try:
            return list(map(float, val_list))
        except:
            return val_list
# -------------------------------------------------------------------- #

def round_up_to_even(f):
    return np.ceil(f / 2.) * 2

# ==============================================#
def generate_RBM_inputs(config, folder): 
    cfg = read_config(config)  # Read config file

    # [INPUT]
    fdr_ds = xr.open_dataset(cfg['OPT']['fdr_nc'])
    print('Done reading flow direction file:', cfg['OPT']['fdr_nc'])

    fdr_df = pd.DataFrame({'hru': fdr_ds.seg_id, 'length': fdr_ds.length, 'next_node': fdr_ds.Tosegment}, columns=['hru', 'length', 'next_node'])

    # Output nc file - energy
    output_energy_nc = cfg['INPUT']['output_energy_nc']

    if os.path.exists(cfg['OPT']['network_csv']) == True:
        Ordered_reaches_final = pd.read_csv(cfg['OPT']['network_csv'])   
    else:             
        # ----------------------------------------------#
        #                                              #
        #     find headwater and outlet grid cells     #
        #                                              #
        # ----------------------------------------------#

        fdr_df['Outlet'] = fdr_df['next_node'].isnull().astype(int)
        fdr_df["Tailwaters"] = fdr_df['hru'].isin(fdr_df["next_node"]).astype(int)

        # ----------------------------------------------------#
        #                                                    #
        #           determine calculation order              #
        #           (headwater -> node)                      #
        #                                                    #
        # ----------------------------------------------------#

        outlet_row = fdr_df.loc[fdr_df['next_node'].isnull().astype(int)==1]
        outlet_index = outlet_row.index[0]

        # Determine order of reaches
        Ordered_reaches = pd.DataFrame(columns=['hru', 'length', 'next_node', 'Outlet', 'Tailwaters', 'subbasin'])
        for value, i in enumerate(list(fdr_df[fdr_df["Tailwaters"] == 0].index)):
            Sub_basin = pd.DataFrame(columns=['hru', 'length', 'next_node', 'Outlet', 'Tailwaters'])
            indexNumbers = []


            indexNumbers.append(i)
            j = fdr_df.loc[fdr_df["hru"] == fdr_df.next_node[i]].index
            indexNumbers.append(j[0])

            while j[0] != outlet_index:
                j = fdr_df.loc[fdr_df["hru"] == fdr_df.next_node[j[0]]].index
                indexNumbers.append(j[0])

            Sub_basin = Sub_basin.append(fdr_df.loc[indexNumbers])
            Sub_basin['subbasin'] = value

            Ordered_reaches = Ordered_reaches.append(Sub_basin)

        # Create final ordered list
        Ordered_reaches["hru_above"] = Ordered_reaches.groupby('subbasin').cumcount() # determine length of the chain

        # Drop long chains of hrus
        seg_order = Ordered_reaches.loc[Ordered_reaches['Outlet'] == 1].sort_values('hru_above', ascending = False)['subbasin']

        Ordered_reaches.subbasin = Ordered_reaches.subbasin.astype("category")
        Ordered_reaches.subbasin.cat.set_categories(seg_order, inplace=True)
        Ordered_list = Ordered_reaches.sort_values(["subbasin", "hru_above"])

        Ordered_reaches_final = Ordered_list.drop_duplicates(subset='hru')

        # Add back in the first duplicate row
        def add_row(x):
            if x['Outlet'].max() < 1:
                next_hru = x.iloc[-1]['next_node']
                last_row = Ordered_reaches_final.loc[Ordered_reaches_final["hru"] == next_hru]
                last_row['subbasin'] = x.iloc[-1]['subbasin']
                last_row['hru_above'] = x.iloc[-1]['hru_above'] + 1
                return x.append(last_row)
            return x
        
        Ordered_reaches_final=Ordered_reaches_final.groupby('subbasin').apply(add_row).reset_index(drop=True)

        # Get all groups with more than 1 row
        Ordered_reaches_final = Ordered_reaches_final.groupby(['subbasin']).filter(lambda x: len(x) >= 2) 

        outlet_subbasin = Ordered_reaches_final.loc[Ordered_reaches_final['Outlet']==1].reset_index(drop=True)
        outlet_subbasin = outlet_subbasin.iloc[0]['subbasin']

        Outlet_basin = Ordered_reaches_final[Ordered_reaches_final['subbasin'] == outlet_subbasin]
        other_basins = Ordered_reaches_final[Ordered_reaches_final['subbasin'] != outlet_subbasin].groupby('subbasin').filter(lambda g: g.isnull().sum().sum() < 1)

        Ordered_reaches_final = pd.concat([other_basins, Outlet_basin])

        Ordered_reaches_final = Ordered_reaches_final.reset_index(drop=True)

        Ordered_reaches_final['node'] = Ordered_reaches_final.index.values + 1

        # GET RIVER LENGTH
        Ordered_reaches_final['seg_length'] = Ordered_reaches_final['length'].groupby(Ordered_reaches_final['subbasin']).transform('sum')
        Ordered_reaches_final['seg_length'] = Ordered_reaches_final['seg_length']*0.000621371 # convert to miles

        Ordered_reaches_final["Upstream_mile"] = Ordered_reaches_final['length'].groupby(Ordered_reaches_final['subbasin']).cumsum() 
        Ordered_reaches_final["Upstream_mile"] = Ordered_reaches_final["Upstream_mile"]*0.000621371 # convert to miles

        # Get river mile
        Ordered_reaches_final['river_mile'] = Ordered_reaches_final['seg_length'] - Ordered_reaches_final['Upstream_mile']

        # Create flow and heat cells
        nreach = len(pd.unique(Ordered_reaches_final['subbasin']))
        heatcell_tot = len(Ordered_reaches_final['hru'])
        flowcell_tot = heatcell_tot - nreach + 1

        # Determine the calculation order
        next_seg_list = []
        for reach in Ordered_reaches_final.subbasin.unique():

            if Ordered_reaches_final.subbasin.unique()[-1] == reach:
                outlet_reach = 0
                to_seg = len(Ordered_reaches_final['hru'])
            else:
                temp_df = Ordered_reaches_final.loc[Ordered_reaches_final.subbasin==reach]
                temp_index = temp_df.index[-1]
                temp=Ordered_reaches_final.loc[Ordered_reaches_final["hru"] == Ordered_reaches_final.next_node[temp_index]]

                if len(temp) > 1:
                    to_seg = int(temp.loc[(temp.river_mile!=0)].node)
                else:
                    to_seg = int(temp.node)

            next_seg_list.append(to_seg)

        dictOfWords = dict(zip(Ordered_reaches_final.subbasin.unique(), next_seg_list))

        Ordered_reaches_final["to_seg"] = Ordered_reaches_final['subbasin'].map(dictOfWords)

        Ordered_reaches_final = Ordered_reaches_final.sort_values(['to_seg','node']).reset_index(drop=True)
        Ordered_reaches_final['node_new'] = Ordered_reaches_final.index.values + 1
        dictOfWords = dict(zip(Ordered_reaches_final.node.values, Ordered_reaches_final.node_new.values))
        Ordered_reaches_final["to_new_seg"] = Ordered_reaches_final['to_seg'].map(dictOfWords)

        # Determine the correct order

        while (all(Ordered_reaches_final.to_new_seg.values[i] <= Ordered_reaches_final.to_new_seg.values[i+1] for i in range(len(Ordered_reaches_final.to_new_seg.values) - 1)) == False):

            # rest columns
            Ordered_reaches_final['node'] = Ordered_reaches_final['node_new']
            Ordered_reaches_final["to_seg"] = Ordered_reaches_final["to_new_seg"]

            # arrange columns
            Ordered_reaches_final = Ordered_reaches_final.sort_values(['to_seg','node']).reset_index(drop=True)
            Ordered_reaches_final['node_new'] = Ordered_reaches_final.index.values + 1
            dictOfWords = dict(zip(Ordered_reaches_final.node.values, Ordered_reaches_final.node_new.values))
            Ordered_reaches_final["to_new_seg"] = Ordered_reaches_final['to_seg'].map(dictOfWords)

        subbasins = Ordered_reaches_final.subbasin.unique()
        dictOfBasins = dict(zip(subbasins, list(range(1, len(subbasins)+1))))
        Ordered_reaches_final["sub_id"] = Ordered_reaches_final['subbasin'].map(dictOfBasins)

        dictOfBasins = dict(zip(list(Ordered_reaches_final.node.values), list(Ordered_reaches_final.sub_id.values)))
        Ordered_reaches_final["to_sub_id"] = Ordered_reaches_final['to_seg'].map(dictOfBasins)       

        Ordered_reaches_final.to_csv(cfg['OPT']['network_csv'])
        
        print('Determine the calculation order of each nodes!')

    # ----------------------------------------------------#
    #                                                    #
    #                    Write output                    #
    #                                                    #
    # ----------------------------------------------------#
       
    # Create flow and heat cells
    nreach = len(pd.unique(Ordered_reaches_final['subbasin']))
    heatcell_tot = len(Ordered_reaches_final['hru'])
    flowcell_tot = heatcell_tot - nreach + 1
    
    input_id = cfg['OPT']['input_id']

    print('Start writing output...')

    site_loc = cfg['OPT']['site_loc']
    data = pd.read_csv(cfg['OPT']['optimization_param'])
    data = data[data['Location'] == site_loc]

    timelag = 0.1

    alpha = data[data.iloc[:,0] == folder].iloc[:,1].values[0]
    beta = data[data.iloc[:,0] == folder].iloc[:,2].values[0]
    gamma = data[data.iloc[:,0] == folder].iloc[:,4].values[0]
    mu = data[data.iloc[:,0] == folder].iloc[:,3].values[0]
    min_velocity = data[data.iloc[:,0] == folder].iloc[:,9].values[0]
    
    # Calculate number of cells necessary
    min_seg = cfg['PARAM']['seg_minflow']
    Ordered_reaches_final['seg_num'] = round_up_to_even(Ordered_reaches_final['length']*0.000621371/(min_velocity/1.5*24)) + min_seg 
    
    # Specify path
    file_path = '/glade/scratch/dblaskey/RBM/RBM_Input/%s/%s_Network'%(folder, input_id)

    # Check whether the specified
    # path exists or not
    if os.path.exists(file_path) == True:
        print("Skip network file for ", input_id, " in ", folder)            
    else:                                    
        print("Creating network file for", input_id, " in ", folder) 

        flow_loc = '/glade/scratch/dblaskey/RBM/RBM_Input/%s/%s_flow'%(folder, input_id)

        network_file = open('/glade/scratch/dblaskey/RBM/RBM_Input/%s/%s_Network'%(folder, input_id), "w")
        network_file.write('RBM parameter preparation network file \n')
        network_file.write('%s\n' % flow_loc)
        network_file.write('%s\n' % cfg['OPT']['rbm_energy_file'])
        network_file.write('/glade/scratch/dblaskey/RBM/routing_dummy.txt \n')
        network_file.write('/glade/scratch/dblaskey/RBM/routing_dummy \n')
        network_file.write('/glade/scratch/dblaskey/RBM/routing_dummy_csv.csv \n')
        network_file.write("%s %s\n" % (cfg['PARAM']['start_date'], cfg['PARAM']['end_date']))

        source = 'False'
        whether_res = 'False'
        Whether_withdraw = 'False'
        res_num = 0
        sub_daily_time_step = cfg['PARAM']['sd_time_step']

        width = (10,10,10,10,10,10,10,10,2)
        items = [nreach, flowcell_tot, heatcell_tot, sub_daily_time_step, source, whether_res, res_num, Whether_withdraw,'\n']
        network_file.write("".join("%*s" % i for i in zip(width, items)))

        width3 = (4,6,4,6,7,6,5,9,5,11,7,10,5,5,2)
        width1 = (5,11,4,9,6,12,8,15,10,2)
        width2 = (8,9,9,9,9,2)

        for value, reach in enumerate(Ordered_reaches_final.subbasin.unique()):
            Subbasin_reaches =  Ordered_reaches_final[Ordered_reaches_final['subbasin']==reach]

            if Ordered_reaches_final.subbasin.unique()[-1] == reach:
                outlet_reach = 0
            else:
                outlet_reach = Subbasin_reaches.to_sub_id.values[-1]

            to_seg = int(Subbasin_reaches.to_new_seg.values[-1] - 1)

            items1 = [len(Subbasin_reaches['hru']), 'Headwaters', int(value + 1), 'TribCell', to_seg, 'Headwaters', int(outlet_reach), 'R.M. =', int(Subbasin_reaches['seg_length'].iloc[-1]), '\n']
            items2 = ['%.2f' % alpha, '%.2f' % beta, '%.4f' % gamma, '%.2f' % mu, '%.4f' % timelag, '\n']
            network_file.write("".join("%*s" % i for i in zip(width1, items1)))
            network_file.write("".join("%*s" % i for i in zip(width2, items2)))

            for i in Subbasin_reaches.index.values:
                rm_node = Subbasin_reaches['river_mile'][i]
                items3 = ['Node', int(Subbasin_reaches['node_new'][i]), 'Row', int(1), 'Column', int(1), 'Lat', int(Ordered_reaches_final['hru'][i]), 'Long', Ordered_reaches_final['next_node'][i], 'R.M. =', '%.2f' % rm_node, int(Ordered_reaches_final['seg_num'][i]), int(0), '\n']
                network_file.write("".join("%*s" % i for i in zip(width3, items3)))
        network_file.close()

        print('Network file:', '%s' % folder)

    Node_list = Ordered_reaches_final[['hru','node_new','subbasin']].astype(int)
    #====================================================#
    # Load and process mizuRoute output data - flow
    #====================================================#
    years = range(int(str(cfg['PARAM']['start_date'])[:4]), int(str(cfg['PARAM']['end_date'])[:4]) + 1)
    
    #=== Read Parameters ===#
    a_d = data[data.iloc[:,0] == folder].iloc[:,5].values[0]
    b_d = data[data.iloc[:,0] == folder].iloc[:,6].values[0]
    a_w = data[data.iloc[:,0] == folder].iloc[:,7].values[0]
    b_w = data[data.iloc[:,0] == folder].iloc[:,8].values[0]
    min_velocity = data[data.iloc[:,0] == folder].iloc[:,9].values[0]
    min_d = data[data.iloc[:,0] == folder].iloc[:,10].values[0]
    min_w = data[data.iloc[:,0] == folder].iloc[:,11].values[0]
    
    #=== Read in Manning's equation data ===#
    #Node_widths = pd.read_csv('/glade/u/home/dblaskey/RBM/manning_param.csv') 
    #Node_widths[Node_widths['Slope']<0.0001] = 0.0001

    for year in years:
        
        # Specify path
        file_path = '/glade/scratch/dblaskey/RBM/RBM_Input/%s/%s_flow_%s'%(folder, input_id, year)

        # Check whether the specified
        # path exists or not
        if os.path.exists(file_path) == True:
            print("Skip flow flie for ", input_id, " in ", folder)            
        else:                                    
            print('Loading and processing mizuRoute output flow data for ' + str(year))

            #=== Load data ===#
            ds_flow = xr.open_dataset(glob.glob(cfg['OPT']['flow_output_nc']+str(year)+'-*')[0])
            da_flow = ds_flow['IRFroutedRunoff']
            da_flow = da_flow.resample(time='D').mean(dim='time')

            da_seg = pd.DataFrame(ds_flow['reachID'], columns = ["hru"])
            da_seg["MR_Output_Order"] = da_seg.index.values 

            Node_map = da_seg.merge(Node_list)
            #Node_widths = Node_map.merge(Node_widths).drop_duplicates()
            
            #=== Convert units ===#
            da_flow = da_flow / 0.028316847 ### cms to cfs

            #=== Set minimum discharge ===#
            da_flow.values[da_flow.values<cfg['RBM_OPTIONS']['min_flow']] = cfg['RBM_OPTIONS']['min_flow']     
            
            flow_loc = '/glade/scratch/dblaskey/RBM/RBM_Input/%s/%s_flow'%(folder, input_id) 

            #=== Calculate flow depth, width and velocity ===#
            da_depth = a_d * pow(da_flow, b_d)  # flow depth [ft]
            da_width = a_w * pow(da_flow, b_w)  # flow width [ft]

            #=== Set minimum width, depth ===# 
            da_width.values[da_width.values< min_w] = min_w
            da_depth.values[da_depth.values< min_d] = min_d
                
            #=== calculate velocity ===# 
            da_velocity = da_flow / da_depth / da_width  # flow velocoty [ft/s]
            
            #=== Set minimum velocity ===# 
            da_velocity.values[da_velocity.values<min_velocity] = min_velocity

            #==== Drop last HRU in every reach ===============#
            Node_map = Node_map.sort_values('node_new')
            final_row = Node_map.iloc[-1] # pull final row to add it back on later
            flow_list = Node_map.groupby('subbasin').apply(lambda x: x.iloc[:-1]).reset_index(drop=True) # drop last HRU
            flow_list = flow_list.append(final_row, ignore_index = True) # add outlet row
            flow_list = flow_list.sort_values('node_new') # arrange so it is in order
            flow_list[['hru', 'MR_Output_Order', 'node_new', 'subbasin']] = flow_list[['hru', 'MR_Output_Order', 'node_new', 'subbasin']].astype(int) 
            
            #====================================================#
            # Rearrange data - flow
            #====================================================#
                        #=== Put data for each HRU into a df ===#
            list_df = [] 
            for i, seg in enumerate(flow_list.MR_Output_Order.values):
                df = pd.DataFrame(index=da_flow.coords['time'].values) # create df
                df['day'] = range(1, len(df)+1)  # day number (1,2,3,...)
                df['cell'] = int(flow_list.loc[flow_list['MR_Output_Order']== seg]['node_new'])  # seg number
                df['Q_in'] = da_flow.loc[:, seg].values  # inflow discharge [cfs]
                df['Q_out'] = df['Q_in']  # outflow discharge [cfs]
                df['Q_diff'] = 0.0  # lateral flow [cfs]
                df['depth'] = da_depth.loc[:, seg].values  # flow depth [ft]
                df['width'] = da_width.loc[:, seg].values  # flow width [ft]
                df['velocity'] = da_velocity.loc[:, seg].values  # flow velocity [ft/s]
                df['Q_local'] = 0.0
                list_df.append(df)
   
            #=== Combine df of all HRUs together, in a single multiindex df (indice: seg; date) ===#
            df_flow = pd.concat(list_df, keys= flow_list.node_new)
            del list_df
            gc.collect()

            #=== Switch order of indice (to: date; cell), then sort ===#
            df_flow = df_flow.reorder_levels([1,0], axis=0)
            df_flow = df_flow.sort_index()

            #=== Writing data to file ===#
            print('Writing flow data to file...')
            np.savetxt(flow_loc+"_"+str(year), df_flow.values, fmt='%d %d %.1f %.1f %.1f %.1f %.1f %.2f %.1f')

    #=== Load data ===#
    spatial_weights = xr.open_dataset(cfg['OPT']['grid_weight'])
    jj = spatial_weights.j_index - 1
    ii = spatial_weights.i_index - 1
    
    #=== Create Node Map ===#
    ds_flow = xr.open_dataset(glob.glob(cfg['OPT']['flow_output_nc']+str(1990)+'-*')[0])
    da_seg = pd.DataFrame(ds_flow['reachID'], columns = ["hru"])
    da_seg["MR_Output_Order"] = da_seg.index.values 

    Node_map = da_seg.merge(Node_list)
    Node_map = Node_map.sort_values('node_new')
    #==== Create energy files ====#
    for year in years:
        
        # Specify path
        file_path = '/glade/scratch/dblaskey/RBM/RBM_Input/%s_energy_%s'%(input_id, year)

        # path exists or not
        if os.path.exists(file_path) == False:                                   
            #====================================================#
            # Load and process RASM/CTSM output data - energy
            #====================================================#
            print('Loading and processing RASM-CTSM output energy data for ' + str(year))

            # Enable dataframes
            Tair = pd.DataFrame(index = np.unique(spatial_weights.IDmask))
            Tair = Tair.rename_axis('COMID')
            vpd = pd.DataFrame(index = np.unique(spatial_weights.IDmask))
            vpd = vpd.rename_axis('COMID')
            Shortwave = pd.DataFrame(index = np.unique(spatial_weights.IDmask))
            Shortwave = Shortwave.rename_axis('COMID')
            Longwave = pd.DataFrame(index = np.unique(spatial_weights.IDmask))
            Longwave = Longwave.rename_axis('COMID')
            Pressure = pd.DataFrame(index = np.unique(spatial_weights.IDmask))
            Pressure = Pressure.rename_axis('COMID')
            Wind = pd.DataFrame(index = np.unique(spatial_weights.IDmask))
            Wind = Wind.rename_axis('COMID')

            if year == 2021:
                energy_ds = xr.open_dataset(cfg['INPUT']['output_energy_nc']+str(year)+'py.nc')
            else:
                energy_ds = xr.open_dataset(cfg['INPUT']['output_energy_nc']+str(year)+'.nc')

            energy_ds = energy_ds.isel(lat = jj,lon = ii)
            energy_ds = energy_ds.resample(time='D').mean(dim='time') 

            # determine the weighted values for each grid cell
            Tair = (spatial_weights.weight*energy_ds.TSA).to_pandas()
            vpd = (spatial_weights.weight*energy_ds.VPD2M).to_pandas()
            Shortwave = (spatial_weights.weight*energy_ds.FSDS).to_pandas()
            Longwave = (spatial_weights.weight*energy_ds.FLDS).to_pandas()
            Pressure = (spatial_weights.weight*energy_ds.PBOT).to_pandas()
            Wind = (spatial_weights.weight*energy_ds.U10).to_pandas()

            # Sum each of the weights of each grid cell the HRU crosses
            Tair["COMID"] = spatial_weights.IDmask
            Tair=Tair.groupby("COMID").sum()
            vpd["COMID"] = spatial_weights.IDmask
            vpd = vpd.groupby("COMID").sum()
            Shortwave["COMID"] = spatial_weights.IDmask
            Shortwave=Shortwave.groupby("COMID").sum()
            Longwave["COMID"] = spatial_weights.IDmask
            Longwave=Longwave.groupby("COMID").sum()
            Pressure["COMID"] = spatial_weights.IDmask
            Pressure=Pressure.groupby("COMID").sum()
            Wind["COMID"] = spatial_weights.IDmask
            Wind=Wind.groupby("COMID").sum()

            print('Done with weights for ' + str(year))

            #=== Set Minimum Wind Speed ===#
            Wind.values[Wind.values< 0.5] = 0.5

            #=== Calculate Vapor Pressure ===#
            # the Clausius-Clapeyron relation: vp = es - vpd = es(T0) X exp(L/Rw*((1/T0) - (1/T))) - vpd (units in mb)
            vp =(6.11 * np.exp(2500000/461.52 * ((1/273.15)-(1/Tair)))) - (vpd/100) 

            #=== Converting units ===#
            Tair = Tair - 273.15 # convert [K] to [C]
            Shortwave = Shortwave * 2.388 * pow(10, -4) # convert [W/m2] to [mm*K/s]
            Longwave = Longwave * 2.388 * pow(10, -4) # convert [W/m2] to [mm*K/s]
            Pressure = Pressure / 100.0 # convert [Pa] to [mb]
            print('Done with unit conversions')

            # Organize for the routing file
            list_df = []
            for i, seg in enumerate(Node_map.MR_Output_Order.values):
                df = pd.DataFrame(index=energy_ds.time.values) # create df
                df['cell'] = int(i+1)  # seg number
                df['Tair'] = Tair.iloc[seg].T.values # Tair
                df['vp'] = vp.iloc[seg].T.values# vp
                df['Shortwave'] = Shortwave.iloc[seg].T.values # Shortwave
                df['Longwave'] = Longwave.iloc[seg].T.values # Longwave
                df['Density'] = 1.225 # Density
                df['Pressure'] = Pressure.iloc[seg].T.values # Pressure
                df['Wind'] = Wind.iloc[seg].T.values # Wind
                list_df.append(df)  

            #=== Combine df of all grid cells together, in a single multiindex df (indice: cell; date) ===#
            df_energy = pd.concat(list_df, keys=Node_map.node_new)

            #=== Switch order of indice (to: date; cell), then sort ===#
            df_energy = df_energy.reorder_levels([1,0], axis=0)
            df_energy = df_energy.sort_index()

            #=== Writing data to file ===#
            print('Writing energy data to file...')
            np.savetxt(cfg['OPT']['rbm_energy_file']+"_"+str(year), df_energy.values, fmt='%d %.1f %.1f %.4f %.4f %.3f %.1f %.1f')

            #=== Clean up memory ===#
            del df_energy
            del list_df
            print('Done with ' + str(year))

print('DONE!')
