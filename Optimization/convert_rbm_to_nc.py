#!/usr/local/anaconda/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sys
import subprocess

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
#==============================================#

def rbm_to_nc(config):
    #==========================================================#
    # read in configuration files
    #==========================================================#
    cfg = read_config(config)
    print(cfg)
    # number of segments
    with open(cfg['INPUT']['spat_file'], 'r') as fp:
        seg_number = len(fp.readlines())

    # start date and end date
    # start_date_str = '{:4d}-{:02d}-{:02d}'.format(cfg['INPUT']['start_date_rbm'][0], \
    #                         cfg['INPUT']['start_date_rbm'][1], \
    #                         cfg['INPUT']['start_date_rbm'][2])
    # end_date_str = '{:4d}-{:02d}-{:02d}'.format(cfg['INPUT']['end_date_rbm'][0], \
    #                         cfg['INPUT']['end_date_rbm'][1], \
    #                         cfg['INPUT']['end_date_rbm'][2])

    # timeperiod = pd.date_range(start_date_str,end_date_str)

    sim_year = int(cfg['PARAM']['sim_year'])

    if sim_year == 2021:
        timeperiod = pd.Series(pd.date_range(start=f'%s-01-01'%sim_year, end=f'%s-09-29'%sim_year, freq='1D'))
    else:
        timeperiod = pd.Series(pd.date_range(start=f'%s-01-01'%sim_year, end=f'%s-12-31'%sim_year, freq='1D'))

    #===========================================================#
    # Select grid cells
    #===========================================================#
    print('Open .Spat file: ', cfg['INPUT']['spat_file'])
    f = open(cfg['INPUT']['spat_file'], 'r')

    list_reach_no = [] # reach number
    list_cell_no = [] # cell number
    list_hru = [] # the hru
    list_segment_no = [] # 1 or 2
    list_execute = [] # 1: will include in the nc file; 0: will not include in the nc file

    for i in range(seg_number):
        line = f.readline().rstrip("\n")
        print(line)
        hru_index = [float(line.split()[4])]
        seg_index = int(float(line.split()[6]))
        reach_index = int(line.split()[0])
        if seg_index == 1:
            if hru_index in list_hru:
                temp_index = [i for i,x in enumerate(list_hru) if x==hru_index]
                reach_max = list_reach_no[np.max(temp_index)]
                if reach_index > reach_max:
                    list_execute.append(1)
                    for j in temp_index:
                        list_execute[j] = 0
                else:
                    list_execute.append(0)
            else:
                list_execute.append(1)
        if seg_index == 2:
            if list_execute[-1]==0:
                list_execute.append(0)
            else:
                list_execute.append(1)
        list_reach_no.append(int(line.split()[0]))
        list_cell_no.append(int(line.split()[1]))
        list_hru.append(hru_index)
        list_segment_no.append(seg_index)

    #===========================================================#
    # Open data sets and define variable range
    #===========================================================#
    print('Open .Temp file ',cfg['INPUT']['temp_file'])
    flowdata = open(cfg['INPUT']['temp_file'], 'r')

    coor_hru = []
    coor_seg = []
    coor_reach = []

    for i in range(seg_number):
        if list_execute[i] == 1:
            coor_hru.append(list_hru[i][0])
            coor_seg.append(list_segment_no[i])

    hru_dim = np.unique(coor_hru)
    seg_dim = np.unique(coor_seg)

    #===========================================================#
    # find index for difference coordinates
    #===========================================================#
    print('Generating index for different variables...')
    index_hru = []
    index_seg = []

    for i in range(seg_number):
        index_hru.append(np.where(hru_dim==list_hru[i])[0][0])
        index_seg.append(np.where(seg_dim==list_segment_no[i])[0][0])

    #===========================================================#
    # Read in file and create dataset
    #===========================================================#
    print('Create empty data data array for each variable...')

    T_water = np.full([len(timeperiod), len(seg_dim), len(hru_dim)],np.NaN)
    T_air = np.full([len(timeperiod), len(seg_dim), len(hru_dim)],np.NaN)
    T_headwater = np.full([len(timeperiod), len(seg_dim), len(hru_dim)],np.NaN)

    print('Loop through every day...')
    for date_index in range(len(timeperiod)):
        print(date_index)
        for i in range(seg_number):
            line = flowdata.readline().rstrip("\n")
            if list_execute[i] == 1:
                T_water[date_index, index_seg[i], index_hru[i]] = float(line.split()[5])
                T_headwater[date_index, index_seg[i], index_hru[i]] = float(line.split()[6])
                T_air[date_index, index_seg[i], index_hru[i]] = float(line.split()[7])

    T_water_da = xr.DataArray(T_water,coords=[timeperiod,seg_dim,hru_dim], dims=['time', 'no_seg', 'hru'])
    T_air_da = xr.DataArray(T_air, coords=[timeperiod,seg_dim,hru_dim], dims=['time', 'no_seg', 'hru'])
    T_head_da = xr.DataArray(T_headwater,coords=[timeperiod,seg_dim,hru_dim], dims=['time', 'no_seg', 'hru'])
    data_sum = xr.merge([T_water_da.rename('T_stream'), T_air_da.rename('T_air'), T_head_da.rename('T_headwater')])    

    print('Write output file to: ', cfg['OUTPUT']['rbm_nc_file'] )
    data_sum.to_netcdf(cfg['OUTPUT']['rbm_nc_file'],format='NETCDF4_CLASSIC')

    # subprocess.call("scp " + cfg['OUTPUT']['rbm_nc_file'] + "", shell=True)
