"""
#################################################################
File: binning.py

Author: Michael Hemming 
Date: 20/03/2024
Description: functions used for binning mooring data
#################################################################
"""

# %% -----------------------------------------------------------------------------------------------
# Determine which computer this script is on

import os
if 'z3526971' in os.getcwd():
    account = 'C:\\Users\\z3526971\\'
else:
    account = 'C:\\Users\\mphem\\'


# %% -------------------------------------------------------------------------
# Import necessary libraries

import numpy as np
import pandas as pd
from scipy.stats import mode
import os
import xarray as xr
from scipy.stats import binned_statistic

# change to scripts path for EAC package loading
os.chdir(account + '\\OneDrive - UNSW\\Work\\EACMooringsVariability\\EACMooringsVariability\\Scripts\\')    
# import Utilities.Functions.mooring as mooring
# import Utilities.Functions.load.directories as options

# %% -------------------------------------------------------------------------
# Functions

def bin_and_fill_data(V, D, t, depth_interval=1, depth_max = 145, time_interval='W', fill_depth=True):
    """
    
    Parameters
    ----------
    V : numpy ndarray
        Variable (e.g. TEMP, U, V, PSAL)
    D : numpy ndarray
        Depth (metres)
    t : numpy datetime64
        time 
    depth_interval : integer
        depth interval for interpolation. The default is 1.
    depth_max : integer, 
        maximum depth for interpolation, but will be 
        affected by deepest measurement. The default is 145.
    time_interval : string
        string argument for pandas date range. The default is 'W'.
    fill_depth : boolean, 
        True of False to interpolate over depth. The default is True.

    Returns
    -------
    binned_ds : xarray dataset
        Dataset containing the binned (no interp) data.
    gridded_ds : xarray dataset
        Dataset containing the gridded (interp) data.

    """
    
    t_finite = t[np.isfinite(t)]
    # get weekly timestamps over record
    weekly_timestamps = pd.date_range(start=np.datetime64(t_finite.min(), 'D'), end=np.datetime64(t_finite.max(), 'D'), freq='W').to_numpy()
    # data containers
    Vbin = []
    Dbin = []
    tbin = []
    Depth_index_range = np.arange(1,depth_max+1,depth_interval)
    Dint = np.ones((len(Depth_index_range),len(weekly_timestamps)))*np.nan
    Vint = np.ones((len(Depth_index_range),len(weekly_timestamps)))*np.nan
    Tint = np.empty((len(Depth_index_range), len(weekly_timestamps)), dtype='datetime64[D]')
    Tint.fill(np.datetime64('NaT'))
    wtn = 0
    for wt in weekly_timestamps:
        # print(wt)
        # get data for the week
        c = np.logical_and(t > wt-np.timedelta64(84,'h'),
                           t <= wt+np.timedelta64(84,'h'))
        Vweek = V[c]
        Dweek = D[c]
        tweek = t[c]
        # remove NaNs
        cnan = np.logical_or(np.isnan(Vweek),np.isnan(Dweek));
        Vweek = Vweek[~cnan]
        Dweek = Dweek[~cnan]
        tweek = tweek[~cnan]
        if len(Vweek) != 0: # Find peaks in data over depth
            count, edges = np.histogram(Dweek,round(max(Dweek)))
            centres = edges[0:-1]+np.round(np.diff(edges),2)[0]
            check_count = count[count != 0]; 
            fcount = count >= np.nanmedian(check_count)
            DepthsWeek = np.round(centres[fcount])
            # check whether some bins too close to one another, remove unnecesary bins
            DepthsWeekFiltered = []
            # for temperature data, select appropriate bins
            if np.nanmean(Vweek)> 10:
                for n in np.arange(1,len(DepthsWeek),1):
                    if DepthsWeek[n] - DepthsWeek[n-1]  > 2:
                        DepthsWeekFiltered.append(DepthsWeek[n-1])
                    if len(DepthsWeekFiltered) != 0 and np.abs(DepthsWeekFiltered[-1] - DepthsWeek[-1]) > 2 and n == len(DepthsWeek)-1:
                        DepthsWeekFiltered.append(DepthsWeek[n])
                DepthsWeekFiltered = np.array(DepthsWeekFiltered)-1
            else:
                DepthsWeekFiltered = DepthsWeek;
            # if only one depth available, specify that depth using the mode depth
            if len(DepthsWeekFiltered) == 0:
                DepthsWeekFiltered = np.array(mode(DepthsWeek)[0])
            # bin data for each recognised depth
            Vbinned = []
            if np.size(DepthsWeekFiltered) > 1:
                for nD in range(len(DepthsWeekFiltered)):
                    cD = np.logical_and(Dweek >= DepthsWeekFiltered[nD]-2,
                                        Dweek < DepthsWeekFiltered[nD]+2);
                    Vbinned.append(np.nanmedian(Vweek[cD]))
            else:
                cD = np.logical_and(Dweek >= DepthsWeekFiltered-2,
                                    Dweek < DepthsWeekFiltered+2);
                Vbinned.append(np.nanmedian(Vweek[cD]))
            Vbin.append(np.array(Vbinned))
            Dbin.append(DepthsWeekFiltered)
            if np.size(DepthsWeekFiltered) > 1:
                tbin.append(np.repeat(wt,len(DepthsWeekFiltered)))
            else:
                tbin.append(wt)
            
            if len(Vbinned) > 1:
                if fill_depth:

                    # interpolate binned profile
                    Dint_ = Depth_index_range
                    Tint_ = np.repeat(wt,len(Dint_))
                    Vint_ = np.interp(Dint_,DepthsWeekFiltered,Vbinned)
                    # remove interpolated data where no data
                    Vint_[Dint_ < min(DepthsWeekFiltered)] = np.nan
                    Vint_[Dint_ > max(DepthsWeekFiltered)] = np.nan
                    # remove interpolated data outside depth threshold
                    Depths2Use = mooring.tools.generate_depth_range(Dbin[-1])                  
                    Vint_[np.where(~np.isin(Dint_, Depths2Use))[0]] = np.nan
                else:
                    Dint_ = np.ones((1,depth_max))*np.nan
                    Vint_ = np.ones((1,depth_max))[0]*np.nan
                    Tint_ = np.repeat(wt,len(Dint_))
            else:
                Dint_ = np.float32(Depth_index_range);
                Vint_ = np.ones((1,depth_max))[0]*np.nan
                Dint_[0:np.int32(Dbin[-1])-1] = np.nan
                Dint_[np.int32(Dbin[-1])::] = np.nan
                Vint_[np.int32(Dbin[-1])-1] = np.array(Vbinned)
                Tint_ = np.repeat(wt,len(Dint_))
                
            Dint[:,wtn] = Dint_
            Vint[:,wtn] = Vint_
            Tint[:,wtn] = Tint_
            wtn = wtn+1
        else:
            # print('No data this week: ' + str(wt))
            wtn = wtn+1
            
        # convert zero-dimensional items in the binned lists for concatenation
        tbin_1d = [np.atleast_1d(item) for item in tbin]
        Dbin_1d = [np.atleast_1d(item) for item in Dbin]
        Vbin_1d = [np.atleast_1d(item) for item in Vbin]
        # concatenate binned data
        tbin_conc = np.concatenate(tbin_1d)
        Dbin_conc = np.concatenate(Dbin_1d)
        Vbin_conc = np.concatenate(Vbin_1d)
        
    
    # create data sets for interpolated and gridded data
    
    # Binned dataset
    # define data with variable attributes
    data_vars = {'VAR':(['TIME'], Vbin_conc), 
                 'DEPTH':(['TIME'], Dbin_conc)}
    # define coordinates
    coords = {'TIME': (['TIME'], tbin_conc)}
    # define global attributes
    attrs = {'comment':'binned data'}
    # create dataset
    binned_ds = xr.Dataset(data_vars=data_vars, 
                    coords=coords, 
                    attrs=attrs)
    
    # gridded dataset
    # define data with variable attributes
    data_vars = {'VAR':(['DEPTH','TIME'], Vint), 
                 'DEPTH_GRID':(['DEPTH','TIME'], Dint)}
    # define coordinates
    coords = {'TIME': (['TIME'], weekly_timestamps),
              'DEPTH': (['DEPTH'], Depth_index_range)}
    # define global attributes
    attrs = {'comment':'gridded data'}
    # create dataset
    gridded_ds = xr.Dataset(data_vars=data_vars, 
                    coords=coords, 
                    attrs=attrs)
    
    return binned_ds, gridded_ds

# %% ---------------------------------------------------------

def bin_data(V, D, t, time_interval='W'):
    """
    Bin a variable and its depth over a specified time interval.

    Parameters
    ----------
    V : numpy ndarray
        variable
    D : numpy ndarray
        depth
    t : numpy ndarray
        Array of timestamps corresponding to the variable and depth data.
    time_interval : str, optional
        Time interval for binning the data. Default is 'W' (weekly).

    Returns
    -------
    tbin : array
        Timestamps for the binned data.
    Dbin : array
        Binned depth values.
    Vbin : array
        Binned variable values.

    """
    
    # get weekly timestamps over record
    weekly_timestamps = pd.date_range(start=np.datetime64(t.min(), 'D'), end=np.datetime64(t.max(), 'D'), freq='W').to_numpy()
    # data containers
    Vbin = []
    Dbin = []
    tbin = []
    for wt in weekly_timestamps:
        # print(wt)
        # get data for the week
        c = np.logical_and(t > wt-np.timedelta64(84,'h'),
                           t <= wt+np.timedelta64(84,'h'))
        Vweek = V[c]
        Dweek = D[c]
        Vbin.append(np.nanmedian(Vweek))
        Dbin.append(np.nanmedian(Dweek))
        tbin.append(wt)
    
    return tbin, Dbin, Vbin

# %% ---------------------------------------------------------

# This one was used for wind - can I use the onw above for wind too?

def binData(time,data,time_bins,method):
    # Convert datetime64 to Unix timestamps for binning function to work
    time = (time - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    time_bins = (time_bins - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')

    # Compute binned statistics (mean) for the data
    binned = pd.DataFrame()
    binned['vals'], _, _ = binned_statistic(time,data,statistic=method, bins=time_bins)
    # Convert the bin edges (time_bins) to datetime64 values
    time_bins_datetime = pd.to_datetime(time_bins, unit='s')
    # Add the timedelta to the datetime values
    binned.index = time_bins_datetime[0:-1] + pd.to_timedelta('12h') 
    
    return binned

# %% ---------------------------------------------------------

def BinGrid_site(varmap,site):

    # check what variables available
    for n in range(len(varmap[site])):
        if 'TEMP' in list(varmap[site][n].variables):
            TEMPind = n
        try:
            if 'VCUR' in list(varmap[site][n].variables):
                VCURind = n
        except:
            VCURind = np.nan 
    #####################################################        
    # bin/grid TEMP
    T = varmap[site][TEMPind].TEMPQCd.values
    D = varmap[site][TEMPind].TEMPQCd_DEPTH.values
    t = varmap[site][TEMPind].TIME.values 
    #########################
    # if surface buoy
    try: 
        T = np.concatenate([T,varmap[site + '_SURFACE'][0].TEMPQCd.values])
        D = np.concatenate([D, varmap[site + '_SURFACE'][0].TEMPQCd_DEPTH.values])
        t = np.concatenate([t, varmap[site + '_SURFACE'][0].TIME.values])
    except KeyError:
        pass
    # get nominal depths
    ND = mooring.getNomDepth(varmap[site][TEMPind], site)
    #########################
    # create data sets
    time_interval='W'
    binnedT, griddedT = mooring.bin_and_fill_data(T, D, t, depth_interval=1, depth_max = 145, time_interval=time_interval, fill_depth=True)
    binnedT.attrs['site'] = site
    binnedT.attrs['sampling'] = time_interval
    binnedT['NOM_DEPTH'] = ND
    griddedT.attrs['site'] = site
    griddedT.attrs['sampling'] = time_interval
    #########################
    # save data sets
    filename = (options.processed_data_path + 
                      site + '_TEMP_binned.nc')
    print('saving: ' + filename)
    binnedT.to_netcdf(filename)
    filename = (options.processed_data_path + 
                      site + '_TEMP_gridded.nc')
    print('saving: ' + filename)
    griddedT.to_netcdf(filename)
    #####################################################        
    # bin/grid VEL
    try:
        V = varmap[site][VCURind].VCURQCd.values
        U = varmap[site][VCURind].UCURQCd.values
        D = varmap[site][VCURind].UVCURQCd_DEPTH.values
        t = varmap[site][VCURind].TIME.values 
    except NameError:
        print('No velocity')
        pass
    #########################
    # create data sets
    try:
        time_interval='W'
        binnedV, griddedV = bin_and_fill_data(V, D, t, depth_interval=1, depth_max = 145, time_interval=time_interval, fill_depth=True)
        binnedU, griddedU = bin_and_fill_data(U, D, t, depth_interval=1, depth_max = 145, time_interval=time_interval, fill_depth=True)
        # combine U and V
        binnedVELV = binnedV; 
        binnedVELV = binnedVELV.rename({'VAR': 'VCUR'})
        binnedVELU = binnedU; 
        binnedVELU = binnedVELU.rename({'VAR': 'UCUR'})
        griddedVEL = griddedV; 
        griddedVEL = griddedVEL.rename({'VAR': 'VCUR'})
        griddedVEL['UCUR'] = griddedU['VAR']    
        # add attributes
        binnedVELV.attrs['site'] = site
        binnedVELV.attrs['sampling'] = time_interval
        binnedVELU.attrs['site'] = site
        binnedVELU.attrs['sampling'] = time_interval        
        griddedVEL.attrs['site'] = site
        griddedVEL.attrs['sampling'] = time_interval
    except NameError:
        binnedVELV = np.nan
        binnedVELU = np.nan
        griddedVEL = np.nan
        pass
    #########################
    # save data sets
    if type(binnedVELV) != float:
        filename = (options.processed_data_path + 
                          site + '_V_binned.nc')
        print('saving: ' + filename)
        binnedVELV.to_netcdf(filename)
        filename = (options.processed_data_path + 
                          site + '_U_binned.nc')
        print('saving: ' + filename)
        binnedVELU.to_netcdf(filename)        
        filename = (options.processed_data_path + 
                          site + '_VEL_gridded.nc')
        print('saving: ' + filename)
        griddedVEL.to_netcdf(filename)
    
    return binnedT,griddedT,binnedVELV,binnedVELU,griddedVEL



