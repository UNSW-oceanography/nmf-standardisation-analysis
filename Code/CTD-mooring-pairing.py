# %% -----------------------------------------------------------------------------------------------
# Ensure working in correct directory (AODN laptop)

import os
repo_path = '/path/to/nmf-analysis-standardisation'
os.chdir(repo_path)

# %% Import packages

import xarray as xr
import numpy as np
import gzip
import glob
import Code.util.holteandtalley.holteandtalley as ht
import Code.util.PickleStuff as ps
import yaml

# %% Site information

with open("sites.yaml", "r") as f:
    s = yaml.safe_load(f)

sites = list(s.keys())
coords = {site: info["coords"] for site, info in s.items()}
depths = {site: info["depths"] for site, info in s.items()}

del s

# %% Load LTSPs

with gzip.open(repo_path + '/Data/TEMP_LTSPs.pkl.gz', 'rb') as f:
    serialized_data = f.read()
LTSPs = ps.pickle.loads(serialized_data)

del serialized_data

# %% Load CTD aggregated file

CTDagg = ps.PickleLoad(repo_path + '/Data/CTDaggData.pickle')

# %% Match Mooring data with CTD profiles

# for each site
# 1. get unique CTD times
# 2. get closest mooring timestamps (rounded to the minute) per CTD profile
# 3. Save the CTD and mooring profile pairs

for s in sites:
    
    print(s)
    
    # get unique CTD times
    un_times = np.unique(CTDagg[s]['TEMP'].TIME.values)
    un_times_mins = un_times.astype('datetime64[m]').astype(np.int64)
    # get unique mooring times
    mooring_times = np.unique(LTSPs[s + '_TEMP']['TEMP_QCd']['TIME'].values)
    mooring_times_mins = mooring_times.astype('datetime64[m]').astype(np.int64)
    
    # extract mooring data
    mooring_T = LTSPs[s + '_TEMP']['TEMP_QCd'].values
    mooring_D = LTSPs[s + '_TEMP']['DEPTH_QCd'].values
    mooring_t = LTSPs[s + '_TEMP']['TEMP_QCd']['TIME'].values

    # create empty lists for looping
    Tnear = []; tnear = []; Dnear = []; nnear = []; 
    TCTD = []; DCTD = []; tCTD = []; nCTD = []; 

    prof_n = 0
    for n in range(len(un_times_mins)):
        # for each unique day with CTD data, get closest mooring data in time (within 15 mins)
        
        ut = un_times_mins[n]
        print(ut)

        # Broadcast and compute absolute differences
        min_diff = np.abs(mooring_times_mins[:, None] - ut)
        # determine where data less than 15 mins away
        if 'VBM100' not in s:
            nearest_c = (min_diff <= 15).squeeze();
        else:
            nearest_c = (min_diff <= 30).squeeze();
        
        if nearest_c.sum() != 0:
            print(ut)
            
            # get nearest mooring data
            nearest_ts = mooring_times[nearest_c];
            nearest_T = mooring_T[np.logical_and(mooring_t >= nearest_ts.min(),
                                                 mooring_t <= nearest_ts.max())]
            nearest_D = mooring_D[np.logical_and(mooring_t >= nearest_ts.min(),
                                                 mooring_t <= nearest_ts.max())]
            nearest_t = mooring_t[np.logical_and(mooring_t >= nearest_ts.min(),
                                                 mooring_t <= nearest_ts.max())]
            
            # get binned mooring profile
            Drange = np.arange(np.nanmin(nearest_D).round(),np.nanmax(nearest_D).round()+2)
            Tbinned = np.ones(np.size(Drange))*np.nan
            Nbinned = np.ones(np.size(Drange))*np.nan
            for nD in Drange:
                c = np.logical_and(nearest_D >= nD-0.5,
                                   nearest_D < nD+0.5)
                Tbinned[int(nD-np.nanmin(nearest_D).round())] = np.nanmedian(nearest_T[c])
                Nbinned[int(nD-np.nanmin(nearest_D).round())] = c.sum()
            c = np.isfinite(Tbinned)
            Tbinned = Tbinned[c]
            Dbinned = Drange[c]
            Nbinned = Nbinned[c]
            # remove binned depths next to each other with less data
            if len(Dbinned) > 2:
                for nD in range(1,len(Dbinned)):
                    try:
                        if Dbinned[nD]-Dbinned[nD-1] < 2:
                            if Nbinned[nD] < Nbinned[nD-1]:
                                Dbinned = np.delete(Dbinned,nD)
                                Tbinned = np.delete(Tbinned,nD)
                                Nbinned = np.delete(Nbinned,nD)
                            else:
                                Dbinned = np.delete(Dbinned,nD-1)
                                Tbinned = np.delete(Tbinned,nD-1)
                                Nbinned = np.delete(Nbinned,nD-1)                                
                    except:
                        pass
            tbinned = np.repeat(np.max(nearest_t),len(Nbinned))
            
            # mooring data
            Tnear.append(Tbinned)
            Dnear.append(Dbinned)
            tnear.append(tbinned)
            nnear.append(np.ones(np.size(tbinned))*prof_n)
            # CTD_data
            CTD_c = CTDagg[s]['TEMP'].TIME.values == un_times[n]
            TCTD.append(CTDagg[s]['TEMP'].values[CTD_c])
            DCTD.append(CTDagg[s]['TEMP']['DEPTH'].values[CTD_c])
            tCTD.append(CTDagg[s]['TEMP']['TIME'].values[CTD_c])
            nCTD.append(np.ones(CTD_c.sum())*prof_n)
            # if 0m temperature NaN, use 1m or 2m
            
            
            # keep n going for next profiles
            prof_n = prof_n + 1
    
    # save information into a dataset
    ds = xr.Dataset({'TEMP_CTD': (['observations_CTD'], np.concatenate(TCTD)),
                    'DEPTH_CTD': (['observations_CTD'], np.concatenate(DCTD)),
                    'TIME_CTD': (['observations_CTD'], np.concatenate(tCTD)),
                    'profile_n_CTD': (['observations_CTD'], np.concatenate(nCTD)),
                    'TEMP_mooring': (['observations_mooring'], np.concatenate(Tnear)),
                    'DEPTH_mooring': (['observations_mooring'], np.concatenate(Dnear)),
                    'TIME_mooring': (['observations_mooring'], np.concatenate(tnear)),
                    'profile_n_mooring': (['observations_mooring'], np.concatenate(nnear))},
                    
                    coords={  'observations_CTD': np.arange(len(np.concatenate(tCTD))),
                        'observations_mooring': np.arange(len(np.concatenate(tnear)))  }    )
    # save dataset
    ds.to_netcdf(repo_path + '/Data' + '/' + s + '_CTD-M_pair.nc')

# %% get nearest SST data for each mooring profile

for s in sites:
    
    T = []; t = []; prof = []; D_prof = []
    
    print(s)
    ds = xr.open_dataset(repo_path + '/Data' + '/' + s + '_CTD-M_pair.nc')
    # for each mooring profile, get timestamp and add satellite data
    for nt in range(len(ds.TIME_mooring.values)):
        date = str(ds.TIME_mooring.values[nt])
        yr = date[0:4]; mn = date[5:7]; dy = date[8:10];
        # create link
        link = ('https://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3S-1d/ngt/' + yr + '/' + yr + mn + dy + '152000-' + 
                'ABOM-L3S_GHRSST-SSTskin-AVHRR_D-1d_night.nc')
        try:
            sat_ds = xr.open_dataset(link)
            # get nearest data_point to mooring (with QC)
            lon = sat_ds.lon.values; lat = sat_ds.lat.values;
            f_lon = np.argmin(abs(lon-coords[s][1]))
            f_lat = np.argmin(abs(lat-coords[s][0]))
            # select temperature
            T.append(sat_ds.sea_surface_temperature[0,f_lat,f_lon].values - sat_ds.sses_bias [0,f_lat,f_lon].values - 273.15 + 0.17)
            # add sat temperature to new data set for concatenation later
            m_prof = ds.profile_n_mooring.values[nt]
            t_prof = ds.TIME_mooring.values[nt]
            t.append(t_prof)
            prof.append(m_prof)
            D_prof.append(0.2)
        except OSError:
            T.append(np.nan)
            m_prof = ds.profile_n_mooring.values[nt]
            t_prof = ds.TIME_mooring.values[nt]
            t.append(t_prof)
            prof.append(m_prof)
            D_prof.append(0.2)
    # save information into a dataset
    ds = xr.Dataset({'TEMP_mooring': (['observations_mooring'], np.array(T)),
                    'DEPTH_mooring': (['observations_mooring'], np.array(D_prof)),
                    'TIME_mooring': (['observations_mooring'], np.array(t)),
                    'profile_n_mooring': (['observations_mooring'], np.array(prof))},
                    
                    coords={'observations_mooring': np.arange(len(np.array(t)))  }    )
  
    # save dataset
    ds.to_netcdf(repo_path + '/Data' + '/' + s + '_CTD-M_pair_satellite.nc')
        
# %% Some additional QC

#==============================================================================
# NRSNSI
ds = xr.open_dataset(repo_path + '/Data/NRSNSI_CTD-M_pair.nc')
# remove ADCP temperatures
for nt in range(np.int32(ds['profile_n_mooring'].max())):
    c = ds['profile_n_mooring'].values == nt
    if np.sum(ds['DEPTH_mooring'].values[c] > 50) > 1:
        depths = ds['DEPTH_mooring'].values[c].copy()
        ADCP_T = np.argwhere(depths == depths.max())
        depths[ADCP_T] = np.nan
        ds['DEPTH_mooring'].values[c] = depths
ds_mod = ds.copy()        
ds.close()
ds_mod.to_netcdf(repo_path + '/Data/NRSNSI_CTD-M_pair_mod.nc')
#==============================================================================
# NRSMAI
ds = xr.open_dataset(repo_path + '/Data/NRSMAI_CTD-M_pair.nc')
# remove ADCP temperatures
for nt in range(np.int32(ds['profile_n_mooring'].max())):
    c = ds['profile_n_mooring'].values == nt
    if np.sum(ds['DEPTH_mooring'].values[c] > 50) > 1:
        depths = ds['DEPTH_mooring'].values[c].copy()
        ADCP_T = np.argwhere(depths == depths.max())
        depths[ADCP_T] = np.nan
        ds['DEPTH_mooring'].values[c] = depths
# remove erroneous profiles
for prof in [8, 14]:
    c = ds['profile_n_mooring'].values == prof
    ds['TEMP_mooring'].values[c] = np.nan
ds_mod = ds.copy()        
ds.close()
ds_mod.to_netcdf(repo_path + '/Data/NRSMAI_CTD-M_pair_mod.nc')
#==============================================================================
# NRSROT
ds = xr.open_dataset(repo_path + '/Data/NRSROT_CTD-M_pair.nc')
# remove earlier configuration (only use profiles with deepest depth > 50m)
c = ds['profile_n_mooring'].values <  45       
ds_mod = ds.copy()
ds_mod['TEMP_mooring'].values[c] = np.nan
ds_mod['DEPTH_mooring'].values[c] = np.nan
# remove erroneous profiles
c = ds_mod['profile_n_mooring'].values == 103
ds_mod['TEMP_mooring'].values[c] = np.nan
c = ds_mod['profile_n_mooring'].values == 60
ds_mod['TEMP_mooring'].values[c] = np.nan
# remove erroneous binned data between 42-47m for specific profiles
for prof in [46, 47, 48, 54, 57, 52, 58, 60, 61, 62, 63, 69, 70, 71, 73, 82, 87, 88, 89, 90, 97, 108, 112, 113, 114, 115, 116, 119, 122, 123, 124, 125, 126, 127, 128, 130, 132, 142, 146, 147, 150, 151]:
    c = np.logical_and(ds['profile_n_mooring'].values == prof,
                       np.logical_and(ds['DEPTH_mooring'].values >= 38,
                                      ds['DEPTH_mooring'].values <= 48)
                       )
    ds_mod['TEMP_mooring'].values[c] = np.nan   
ds.close()    
ds_mod.to_netcdf(repo_path + '/Data/NRSROT_CTD-M_pair_mod.nc')
#==============================================================================
# NRSKAI
ds = xr.open_dataset(repo_path + '/Data/NRSKAI_CTD-M_pair.nc')
# remove earlier configuration (only use profiles with deepest depth > 50m)
c = ds['TIME_mooring'].values < np.datetime64('2011-10-01')
ds_mod = ds.copy()
ds_mod['TEMP_mooring'].values[c] = np.nan
ds_mod['DEPTH_mooring'].values[c] = np.nan
ds_mod.to_netcdf(repo_path + '/Data/NRSKAI_CTD-M_pair_mod.nc')
#==============================================================================
# PH100
ds = xr.open_dataset(repo_path + '/Data/PH100_CTD-M_pair.nc')
# remove erroneous profiles
for prof in [25, 120, 123, 136]:
    c = ds['profile_n_mooring'].values == prof
    ds['TEMP_mooring'].values[c] = np.nan
ds_mod = ds.copy()        
ds.close()
ds_mod.to_netcdf(repo_path + '/Data/PH100_CTD-M_pair_mod.nc')

# %% If 0m CTD data is missing, extrapolate from 1m, 2m etc. to the surface

files = glob.glob(repo_path + '/Data/' + '*_CTD-M_pair.nc')

for file in files:
    print(file)
    if 'NRSROT' in file or 'NRSKAI' in file or 'NRSNSI' in file or 'NRSMAI' in file or 'PH100' in file:
        fileds = file.replace('_CTD-M_pair.nc','_CTD-M_pair_mod.nc')
    else:   
        fileds = file
        
    with xr.open_dataset(fileds) as ds:    
        un_profs = np.int32(np.unique(ds['profile_n_CTD']))
        for up in un_profs:
            c = ds['profile_n_CTD'].copy().values == up
            c0 = np.argwhere(ds['DEPTH_CTD'].copy().values[c] == 0)
            temp_array = np.array(ds['TEMP_CTD'].copy().values)[c].copy()
            depth_array = np.round(ds['DEPTH_CTD'].values[c]).copy()
            # Compute T0, T1, T2 appropriately:
            try:
                T0 = temp_array[ds['DEPTH_CTD'].values[c] == 0].copy()[0]
                T1 = temp_array[ds['DEPTH_CTD'].values[c] == 1].copy()[0]
                T2 = temp_array[ds['DEPTH_CTD'].values[c] == 2].copy()[0]   
                
                if np.isnan(T0):
                    if np.isfinite(T1):
                        temp_array[c0] = np.array(T1)
                    else:
                        temp_array[c0] = T2   
                    
                ds['TEMP_CTD'].data[c] = temp_array  # Update the DataArray with the modified array
                
                   
            except IndexError:
                pass
            
        ds_mod = ds.copy()
    ds_mod.to_netcdf(file.replace('_CTD-M_pair.nc','_CTD-M_pair_QCd.nc'))
