# %% -------------------------------------------------------------
# Packages

import xarray as xr
import numpy as np
import requests
import re
import PickleStuff as ps
import time
from dask import delayed, compute

# %% -------------------------------------------------------------
# download and concatenate CTD profiles

# function to get file names
def get_thredds_filenames(link):
    # scrape thredds server
    res = requests.get(link)
    txt = res.text
    # identify file names
    start_file = []
    end_file = []
    files = []
    # get start and end locations of file names in text
    for m in re.finditer('IMOS_ANMN', txt):
        start_file.append(m.start())
    for m in re.finditer('.nc', txt):
        end_file.append(m.end())  
    # get file names
    for n_file in range(len(start_file)):
        files.append(txt[start_file[n_file]:end_file[n_file]])
    # Only use unique file names
    files = np.unique(files)
    locs = []
    dates = []
    for n in range(len(files)):
        locs.append(link.replace('catalog.html','') +  files[n])
        locs[n] = locs[n].replace('catalog','dodsC')
        f_und = []
        for m in re.finditer('_', files[n]):
            f_und.append(m.end()) 
        fus = f_und[2]
        
        yr = files[n][fus:fus+4]
        mn = files[n][fus+4:fus+6]
        dy = files[n][fus+6:fus+8]
        hr = files[n][fus+9:fus+11]
        mins = files[n][fus+11:fus+13] 
        dates.append(np.datetime64(yr + '-' + mn + '-' + dy + ' ' + hr + ':' + mins))
        
    return files, locs, dates

def load_single_dataset(filepath, variable):
    print(f"Loading: {filepath}")
    ds = xr.open_dataset(filepath)
    var = ds[variable]
    qc = ds[variable + '_quality_control'] if (variable + '_quality_control') in ds.variables else None
    return var, qc

def remove_instance(ds):
    if type(ds) == tuple:
        ds = ds[0]
    if 'INSTANCE' in ds.coords:
        ds = ds.drop("INSTANCE")
        ds = ds.squeeze("INSTANCE")
        
    return ds

def AggregateProfiles(link, variable):
    # Get filenames, locations, and dates from your THREDDS catalog
    files, locs, dates = get_thredds_filenames(link)
    
    # Create delayed objects for each file load
    delayed_loads = [delayed(load_single_dataset)(l, variable) for l in locs]
    
    # Trigger parallel loading
    results = compute(*delayed_loads)
    
    # Unpack the results into two lists
    CTD_data = []
    CTD_data_QC = []
    for var, qc in results:
        CTD_data.append(remove_instance(var))
        if qc is not None:
            CTD_data_QC.append(remove_instance(qc))
    
    print('Aggregating CTD profiles.')
    # Concatenate along the 'DEPTH' dimension
    Profs = xr.concat(CTD_data, dim='DEPTH', coords='minimal', compat='override')
    Profs_QC = xr.concat(CTD_data_QC, dim='DEPTH', coords='minimal', compat='override') if CTD_data_QC else None
    
    # Build a dictionary with the aggregated profiles.
    AggProfs = {variable: Profs}
    if Profs_QC is not None:
        AggProfs[variable + '_quality_control'] = Profs_QC
    
    print('Done.')
    return AggProfs


# %% -------------------------------------------------------------
# Use the code to get concatenated data sets

sites = ['NRSNSI','CH050', 'CH070', 'CH100', 'SYD100', 'SYD140', 'PH100', 'BMP070', 'BMP090', 'BMP120', 'NRSMAI', 'NRSKAI', 'NRSROT', 'NRSYON']
regions = ['NRS','NSW', 'NSW', 'NSW', 'NSW', 'NSW', 'NSW', 'NSW', 'NSW', 'NSW', 'NRS', 'NRS', 'NRS', 'NRS']

CTDagg = {}
for n in range(len(sites)):
    print(sites[n])
    start = time.perf_counter()
    link = f'https://thredds.aodn.org.au/thredds/catalog/IMOS/ANMN/{regions[n]}/{sites[n]}/Biogeochem_profiles/catalog.html'
    
    AggProfs = AggregateProfiles(link,'TEMP')
    CTDagg[sites[n]] = {}
    CTDagg[sites[n]]['TEMP']  = AggProfs['TEMP']
    CTDagg[sites[n]]['TEMP_quality_control']  = AggProfs['TEMP_quality_control']
    
    elapsed = time.perf_counter() - start
    print(f"Iteration {sites[n]} took {elapsed:.2f} seconds")

ps.PickleSave('/home/mphemming/Documents/GitHub/anmn-analysis-standardisation/Data/CTDaggData.pickle',CTDagg)

# %% -------------------------------------------------------------
# Notes

# profiles concatenated by 'DEPTH', not 'OBSERVATION' like LTSPs
# No concatenation of global attributes or instrument information etc
# There's no profile number variable
# Only tried for TEMP and PSAL at PH and MAI








