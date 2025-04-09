# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 14:49:52 2024

@author: mphem
"""

# %% -----------------------------------------------------------------------------------------------
# Ensure working in correct directory (AODN laptop)

import os
repo_path = '/path/to/nmf-analysis-standardisation'
os.chdir(repo_path)

# %% Import packages

import xarray as xr
import os
import requests
import re
import numpy as np
import Code.util.CTD_AGG.aggregated_profiles as agg
import Code.util.PickleStuff as ps

# %% Get data

sites = ['NRSNSI', 'CH050', 'CH070', 'CH100', 'SYD100', 'SYD140', 'PH100', 'BMP070', 'BMP090', 'BMP120', 'NRSMAI', 'NRSKAI', 'NRSROT', 'NRSYON', 'WATR50','VBM100']

CTDdata = {}
for s in sites:

    # get correct link for downloading data
    link = 0
    if 'NRS' in s:
        link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/NRS/' + s + '/Biogeochem_profiles/catalog.html'
    if s == 'WATR50':
        link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/WA/' + s + '/Biogeochem_profiles/catalog.html'
    if link == 0:
        link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/NSW/' + s + '/Biogeochem_profiles/catalog.html'
    if 'VBM' in s:
        link = 'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/SA/' + s + '/Biogeochem_profiles/catalog.html'
    
    print(link)
    
    # get data
    CTDdata[s] = agg.AggregateProfiles(link,'TEMP')
        
        
# %% save data as a pickle

ps.PickleSave(repo_path + '/Data/CTDaggData.pickle',CTDdata)
