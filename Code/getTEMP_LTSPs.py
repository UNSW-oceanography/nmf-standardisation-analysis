# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:14:28 2024

@author: mphemming
"""

# %% -----------------------------------------------------------------------------------------------
# Ensure working in correct directory

import os
repo_path = '/path/to/nmf-analysis-standardisation'
os.chdir(repo_path)

# %% Import packages

import xarray as xr
import os
import requests
import re
import numpy as np
import gzip
import Code.util.PickleStuff as ps
import yaml

# %% Site information

with open("sites.yaml", "r") as f:
    s = yaml.safe_load(f)

sites = list(s.keys())

del s

# %% functions to get LTSPs

# Get filenames
def getLTSPs(sites):
    files = []
    data_links = []
    site_info = []
    for s in sites:
        # get correct link for downloading data
        link = 0
        if 'NRS' in s:
            link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/NRS/' + s + '/aggregated_timeseries/catalog.html'
        if 'NSW' in s:
            link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/NSW/' + s + '/aggregated_timeseries/catalog.html'
        if 'VBM' in s:
            link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/SA/' + s + '/aggregated_timeseries/catalog.html'
        print(link)
        
        # Send a GET request to the URL
        response = requests.get(link)
        # Use regular expression to find the NetCDF file names in the HTML content
        nc_file_names = re.findall(r"<a href='([^']+\.nc)'><tt>[^<]+\.nc</tt></a>", response.text)
        # get file names
        for n in range(len(nc_file_names)):
            f = nc_file_names[n].find('IMOS_ANMN')
            nc_file_names[n] = nc_file_names[n][f::]
            # get OPenDAP links
            if 'NSW'in s:
                data_links.append('https://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NSW/' + s + 
                 '/aggregated_timeseries/' + nc_file_names[n])
            if 'NRS' in s:
                data_links.append('https://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/' + s + 
                 '/aggregated_timeseries/' + nc_file_names[n])   
            if 'VBM' in s:
                data_links.append('https://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/SA/' + s + 
                 '/aggregated_timeseries/' + nc_file_names[n])
            site_info.append(s)
        files.append(nc_file_names)
    return files, data_links

# load LTSPs
def loadLTSPs(data_links,variable):
    data = {}
    for dls in data_links:
        if variable in dls:
            print(dls)
            ds = xr.open_dataset(dls) 
            data[ds.site_code + '_' + list(ds.variables)[0]] = ds.load()
        
    return data

# %% get LTSPs

_, data_links = getLTSPs(sites)
data2 = loadLTSPs(data_links,'TEMP')

# %% Some QC (only use flag == 1)

data_QCd = {}

for s in sites:
    print(s)
    ds = data[s + '_TEMP']
    # get data
    T = ds['TEMP'].values
    TQC = ds['TEMP_quality_control'].values
    D = ds['DEPTH'].values
    DQC = ds['DEPTH_quality_control'].values
    # apply QC
    ds['TEMP_QCd'] = ds['TEMP']
    ds['TEMP'].values[TQC != 1] = np.nan
    ds['DEPTH_QCd'] = ds['DEPTH']
    ds['DEPTH'].values[DQC != 1] = np.nan   
    # remove unneccesary variables
    variables_to_keep = ['TEMP_QCd', 'DEPTH_QCd', 'TIME']
    filtered_ds = {key: value for key, value in ds.items() if key in variables_to_keep}
    # save in new dataset
    data_QCd[s + '_TEMP'] = filtered_ds
    
# %% save data sets as a pickle

# Serialize the data
serialized_data = ps.pickle.dumps(data_QCd)

# Define the filename for the compressed pickle file
filename = repo_path + r'Data/TEMP_LTSPs.pkl.gz'

# Compress and save the serialized data to a file
with gzip.open(filename, 'wb') as f:
    f.write(serialized_data)

# %%
