# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 13:46:07 2024

@author: mphem
"""

# NOTES

"""
1. Get the LTSPs for mooring sites around Australia
2. Get all of the CTD profiles for the chosen sites
3. Extract the matching mooring-CTD profiles (nearest ones in time, but if more than n minutes discard). Create an save new variables as pickle.
4. For each chosen mooring profile, select SST nearest in time/space and add to profile. Consider the time diff between satellite and mooring. 
5. Estimate MLD for CTD and mooring profile
6. Create plots

"""

# %% -----------------------------------------------------------------------------------------------
# Ensure working in correct directory (AODN laptop)

import os
repo_path = '/path/to/nmf-analysis-standardisation'
os.chdir(repo_path)

# %% Import packages

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
import Code.util.holteandtalley.holteandtalley as ht
import Code.util.PickleStuff as ps
import yaml
from sklearn.metrics import r2_score, mean_squared_error

# %% Site information

with open("sites.yaml", "r") as f:
    s = yaml.safe_load(f)

sites = list(s.keys())
coords = {site: info["coords"] for site, info in s.items()}
depths = {site: info["depths"] for site, info in s.items()}

del s

# %% load pair nc files

files = glob.glob(repo_path + '/Data/' + '*_CTD-M_pair_QCd.nc')
files_sat = glob.glob(repo_path + '/Data/' + '*_CTD-M_pair_satellite.nc')

# load and combine data for each site
pairs = []; pairs_sat = []
for file in files:
    pairs.append(xr.open_dataset(file))
    pairs_sat.append(xr.open_dataset(file.replace('_CTD-M_pair_QCd.nc','_CTD-M_pair_satellite.nc')))
    
# combine
pairs_ds = {} 
for n in range(len(files)):
    s = files[n].split('/')[7].split('_')[0]; 
    pairs_ds[s + '_TEMP_mooring'] = np.concatenate([pairs[n]['TEMP_mooring'].values,
                                                    pairs_sat[n]['TEMP_mooring'].values])
    pairs_ds[s + '_profile_n_mooring'] = np.concatenate([pairs[n]['profile_n_mooring'].values,
                                                         pairs_sat[n]['profile_n_mooring'].values])
    pairs_ds[s + '_TIME_mooring'] = np.concatenate([pairs[n]['TIME_mooring'].values,
                                                    pairs_sat[n]['TIME_mooring'].values])
    pairs_ds[s + '_DEPTH_mooring'] = np.concatenate([pairs[n]['DEPTH_mooring'].values,
                                                     pairs_sat[n]['DEPTH_mooring'].values])

# %% Explore data availability of sites

for s in sites:
     
    print(s)
    
    # get mooring data
    Tm = pairs_ds[s + '_TEMP_mooring']
    Dm = pairs_ds[s + '_DEPTH_mooring']     
    pnm = pairs_ds[s + '_profile_n_mooring']
    tm = pairs_ds[s + '_TIME_mooring']    
    un_profs = np.int32(np.unique(pnm))
    # get CTD data
    CTDds = xr.open_dataset(repo_path + '/Data/' + s + '_CTD-M_pair.nc')
    TCTD = CTDds['TEMP_CTD'].values
    DCTD = CTDds['DEPTH_CTD'].values
    pnCTD = CTDds['profile_n_CTD'].values
    tCTD = CTDds['TIME_CTD'].values
    
    ###############################################################
    # Plot Temp vs Depth for all profiles
    # Create a figure with a custom size
    plt.figure(figsize=(8, 12))
    # Plot the temperature vs. depth data using dot markers
    plt.plot(Tm, Dm, '.', color='k', markersize=8)
    # Label the axes and add a title
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Depth (m)")
    plt.title(s + ' - paired profiles')
    # Invert the y-axis so that greater depth values appear lower on the plot
    plt.gca().invert_yaxis()
    # Add gridlines for better readability
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    # Adjust layout to prevent clipping of labels
    plt.tight_layout()
    
    plt.savefig(repo_path + '/Plots/DataAvailability/' + s + 'Mooring_Temp_vs_Depth.png', dpi=300)
    plt.close()

    ###############################################################
    # Plot Time vs Depth for all profiles
    # Create a figure with a custom size
    plt.figure(figsize=(12, 8))
    # Plot the temperature vs. depth data using dot markers
    plt.plot(tm, Dm, '.', color='k', markersize=8)
    # Label the axes and add a title
    plt.ylabel("Depth (m)")
    plt.title(s + ' - paired profiles')
    # Invert the y-axis so that greater depth values appear lower on the plot
    plt.gca().invert_yaxis()
    # Add gridlines for better readability
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    # Adjust layout to prevent clipping of labels
    plt.tight_layout()
    
    plt.savefig(repo_path + '/Plots/DataAvailability/' + s + 'Mooring_Time_vs_Depth.png', dpi=300)
    plt.close()
    
# %% Estimate MLD and Create comparison plots

for s in sites:
     
    print(s)
    
    # get mooring data
    Tm = pairs_ds[s + '_TEMP_mooring']
    Dm = pairs_ds[s + '_DEPTH_mooring']
    pnm = pairs_ds[s + '_profile_n_mooring']
    tm = pairs_ds[s + '_TIME_mooring']
    un_profs = np.int32(np.unique(pnm))
    # get CTD data
    CTDds = xr.open_dataset(repo_path + '/Data/' + s + '_CTD-M_pair.nc')
    TCTD = CTDds['TEMP_CTD'].values
    DCTD = CTDds['DEPTH_CTD'].values
    pnCTD = CTDds['profile_n_CTD'].values
    tCTD = CTDds['TIME_CTD'].values


    for up in un_profs:
        # print(up)
        # select data for profile
        cCTD = pnCTD == up;
        cM = pnm == up;
        DrangeM = np.arange(np.nanmin(Dm[cM]).round(),np.nanmax(Dm[cM]).round())
        c = np.argsort(Dm[cM])
        
        # Get inteprolated Mooring temperature line
        Tline = Tm[cM][c]; Dline = Dm[cM][c];
        if np.isfinite(Tline).sum() > 1:
            Mline = np.interp(DrangeM,Dline[np.isfinite(Tline)],Tline[np.isfinite(Tline)])
            
            if np.isfinite(Tm[cM]).sum() > 1 and np.isfinite(Mline).sum() > 1 and np.nanmin(DrangeM[np.isfinite(Mline)]) < 10:

                # interpolate gaps in vertical between min and max depth
                minD = np.argmin(DrangeM[np.isfinite(Mline)])
                maxD = np.argmax(DrangeM[np.isfinite(Mline)])
                if np.any(np.isnan(Mline[minD:maxD]) > 0):
                    new_DrangeM = np.arange(DrangeM[np.isfinite(Mline)][minD],DrangeM[np.isfinite(Mline)][maxD])
                    Mline = np.interp(new_DrangeM,DrangeM[np.isfinite(Mline)],Mline[np.isfinite(Mline)])
                    DrangeM = new_DrangeM
                
                # estimate MLD
                hCTD = ht.HolteAndTalley(DCTD[cCTD][np.isfinite(TCTD[cCTD])],TCTD[cCTD][np.isfinite(TCTD[cCTD])]).temp.pressures[
                    ht.HolteAndTalley(DCTD[cCTD][np.isfinite(TCTD[cCTD])],TCTD[cCTD][np.isfinite(TCTD[cCTD])]).temp.calculateTTMLD()]
                hM = ht.HolteAndTalley(DrangeM[np.isfinite(Mline)],Mline[np.isfinite(Mline)]).temp.pressures[
                    ht.HolteAndTalley(DrangeM[np.isfinite(Mline)],Mline[np.isfinite(Mline)]).temp.calculateTTMLD()]
                
                # create figure
                plt.figure(figsize=(4, 8))
                plt.plot(TCTD[cCTD],DCTD[cCTD],zorder=0,label='CTD profile')
                plt.plot(Mline,DrangeM,c='black',zorder=0,label='Mooring Interpolated')
                plt.scatter(Tm[cM],Dm[cM],facecolors='white', edgecolors='black',zorder=1,label='Mooring observed')
                plt.scatter(np.interp(hM,DrangeM,Mline),hM,facecolors='black', edgecolors='black',zorder=1,label='Mooring MLD')
                plt.scatter(np.interp(hCTD,DCTD[cCTD],TCTD[cCTD]),hCTD,facecolors='blue', edgecolors='black',zorder=1,label='CTD MLD')
                stitle = s
                stitle = stitle.replace('\\','')
                plt.title(stitle + ' profile ' + str(up))
                plt.gca().invert_yaxis()
                plt.grid()
                plt.legend(loc='upper left')
                plt.xlim([np.nanmin(TCTD[cCTD])-1, np.nanmax(TCTD[cCTD])+1])
        
                plt.annotate('CTD time: ' + str(tCTD[cCTD][0])[0:16], xy=(0.7, 0.9), xycoords='axes fraction', xytext=(0.3, 0.15),textcoords='axes fraction')
                plt.annotate('Mooring time: ' + str(tm[cM][0])[0:16], xy=(0.7, 0.9), xycoords='axes fraction', xytext=(0.2, 0.1),textcoords='axes fraction')
        
                plt.savefig((repo_path + '/Plots/ProfileComparisons/' +
                            s + '/' + s + '_prof_' + str(up) + '.png'),dpi=300)
                plt.close()

# %% get mean errors between interpolated mooring profiles and CTD profiles 

# Define a function to map month to season
def season_from_date(dt):
    """Return the season for a given Python datetime object using Northern Hemisphere definitions."""
    month = []
    seasons = []
    for n in range(len(dt)):
        month.append(dt[n].month)
        if dt[n].month in [12, 1, 2]:
            seasons.append('summer')
        if dt[n].month in [3, 4, 5]:
            seasons.append('autumn')
        if dt[n].month in [6, 7, 8]:
            seasons.append('winter')
        if dt[n].month in [9, 10, 11]:
            seasons.append('spring')
    return np.array(month), np.array(seasons)

def get_errors(pairs_ds, sites, season):

    errors = {}

    for s in sites:
        
        print(s)
        
        # get mooring data
        Tm = pairs_ds[s + '_TEMP_mooring']
        Dm = pairs_ds[s + '_DEPTH_mooring']
        pnm = pairs_ds[s + '_profile_n_mooring']
        tm = pairs_ds[s + '_TIME_mooring']
        _, snsM = season_from_date([pd.Timestamp(d).to_pydatetime() for d in tm])
        un_profs = np.int32(np.unique(pnm))
        # get CTD datasns = season_from_date([pd.Timestamp(d).to_pydatetime() for d in tm])

        CTDds = xr.open_dataset(repo_path + '/Data/' + 
                                s + '_CTD-M_pair_QCd.nc')
        TCTD = CTDds['TEMP_CTD'].values
        DCTD = CTDds['DEPTH_CTD'].values
        pnCTD = CTDds['profile_n_CTD'].values
        tCTD = CTDds['TIME_CTD'].values
        _, snsCTD = season_from_date([pd.Timestamp(d).to_pydatetime() for d in tCTD])
        
        
        # exclude profiles based on season
        if season != 'all':
            # moorings
            c = snsM == season
            Tm = Tm[c]
            Dm = Dm[c]
            pnm = pnm[c]
            tm = tm[c]
            # CTD
            c = snsCTD == season
            TCTD = TCTD[c]
            DCTD = DCTD[c]
            pnCTD = pnCTD[c]
            tCTD = tCTD[c]
            
        
        DrangeM = np.arange(0,140)
        MCTD_diff = np.ones([len(DrangeM),len(un_profs)])
        Mint = np.ones([len(DrangeM),len(un_profs)])
        MDint = np.ones([len(DrangeM),len(un_profs)])
        CTDint = np.ones([len(DrangeM),len(un_profs)])
        CTDDint = np.ones([len(DrangeM),len(un_profs)])
        MCTD_diff_sat = np.ones([len(DrangeM),len(un_profs)])
        
        nprofs_used = 0
        nprofs_used_sat = 0
        
        
        for n in range(len(un_profs)):
            up = un_profs[n]
            # print(up)
            # select data for profile
            cCTD = pnCTD == up;
            cM = pnm == up;
            c = np.argsort(Dm[cM])
            # select mooring and CTD data
            MD = Dm[cM][c]; MT = Tm[cM][c]
            CTDD = DCTD[cCTD]; CTDT = TCTD[cCTD]
            satT = MT[MD == 0.2]
            # checks
            sat_check = np.any(MD == 0.2)
            sat_ctd_check = []
            depths_check = np.sum(np.isfinite(MD)) >= len(depths[s])
            if sat_check:
                sat_CTD_diff = MT[MD == 0.2][0] - np.nanmedian(CTDT[CTDD <= CTDD.min()+5])
                if np.abs(sat_CTD_diff) < 1:
                    print(np.abs(sat_CTD_diff))
                    sat_ctd_check = True
                else:
                    sat_ctd_check = False
                    print(up)
            CTD_check = np.sum(np.isfinite(CTDT)) >= 30; # at least 30 data points in CTD profile
                    
            # if no satellite data, change sat_ctd_check = 1 to include non-sat data
            if type(sat_ctd_check) == list:
                sat_ctd_check = True        
                    
            # All profiles except those where sat data is > 1 deg C diff
            if depths_check and CTD_check and sat_ctd_check:
                nprofs_used += 1
                # Get inteprolated Mooring temperature line
                Mline = np.interp(DrangeM,MD,MT)
                Mline[DrangeM > np.nanmax(MD)] = np.nan
                Mint[:,n] = Mline
                MDint[:,n] = DrangeM
                # get CTD profile over same range
                CTDline = np.interp(DrangeM,DCTD[cCTD],TCTD[cCTD])
                CTDline[DrangeM > np.nanmax(DCTD[cCTD])] = np.nan
                CTDint[:,n] = CTDline
                CTDDint[:,n] = DrangeM
                # get diff and save in file
                MCTD_diff[:,n] = np.abs(Mline-CTDline)
            else:
                MCTD_diff[:,n] = np.nan                
                        
            # Only satellite profiles (when diff < 1 deg C)
            if depths_check and sat_check and sat_ctd_check and CTD_check:
                nprofs_used_sat += 1
                # Get inteprolated Mooring temperature line
                Mline = np.interp(DrangeM,MD,MT)
                Mline[DrangeM > np.nanmax(MD)] = np.nan
                # get CTD profile over dame range
                CTDline = np.interp(DrangeM,DCTD[cCTD],TCTD[cCTD])
                CTDline[DrangeM > np.nanmax(DCTD[cCTD])] = np.nan
                # get diff and save in file
                MCTD_diff_sat[:,n] = np.abs(Mline-CTDline)
            else:
                MCTD_diff_sat[:,n] = np.nan
                
            del sat_check, sat_ctd_check, CTD_check, depths_check
                
        print(s + ' profiles used: ' + str(nprofs_used))
        print(s + ' profiles used satellite: ' + str(nprofs_used_sat))
        
        MCTD_diff_mean = np.nanmean(MCTD_diff,1)
        MCTD_diff_std = np.nanstd(MCTD_diff,1)  
        MCTD_diff_mean_sat = np.nanmean(MCTD_diff_sat,1)
        MCTD_diff_std_sat = np.nanstd(MCTD_diff_sat,1)  
        
        errors[s + '_mooring_temp'] = Mint
        errors[s + '_mooring_depth'] = MDint
        errors[s + '_CTD_temp'] = CTDint
        errors[s + '_CTD_depth'] = CTDDint
        errors[s + '_diff'] = MCTD_diff
        errors[s + '_diff_mean'] = MCTD_diff_mean
        errors[s + '_diff_std'] = MCTD_diff_std
        errors[s + '_diff_sat'] = MCTD_diff_sat
        errors[s + '_diff_mean_sat'] = MCTD_diff_mean_sat
        errors[s + '_diff_std_sat'] = MCTD_diff_std_sat
        errors[s + '_n_profiles_available'] = un_profs.max()
        errors[s + '_n_profiles_used'] = nprofs_used 
        errors[s + '_n_profiles_used_sat'] = nprofs_used_sat
        
    return errors

# get errors all profiles
errors = get_errors(pairs_ds, sites, 'all')
errors_summer = get_errors(pairs_ds, sites, 'summer')
errors_winter = get_errors(pairs_ds, sites, 'winter')

# %% Create figure for mean error

def plot_error(errors, season):

    plt.figure(figsize=(8, 8))
    
    # Plot the shaded area representing ± standard deviation
    plt.fill_betweenx(np.arange(0,140), errors['NRSNSI_diff_mean'] - errors['NRSNSI_diff_std'], 
                    errors['NRSNSI_diff_mean'] + errors['NRSNSI_diff_std'], 
                    color='blue', alpha=0.1)
    plt.fill_betweenx(np.arange(0,140), errors['NRSMAI_diff_mean'] - errors['NRSMAI_diff_std'], 
                    errors['NRSMAI_diff_mean'] + errors['NRSMAI_diff_std'], 
                    color='orange', alpha=0.1)
    plt.fill_betweenx(np.arange(0,140), errors['NRSROT_diff_mean'] - errors['NRSROT_diff_std'], 
                    errors['NRSROT_diff_mean'] + errors['NRSROT_diff_std'], 
                    color='green', alpha=0.1)
    plt.fill_betweenx(np.arange(0,140), errors['NRSKAI_diff_mean'] - errors['NRSKAI_diff_std'], 
                    errors['NRSKAI_diff_mean'] + errors['NRSKAI_diff_std'], 
                    color='red', alpha=0.1)
    plt.fill_betweenx(np.arange(0,140), errors['PH100_diff_mean'] - errors['PH100_diff_std'], 
                    errors['PH100_diff_mean'] + errors['PH100_diff_std'], 
                    color='purple', alpha=0.1)

    # plot error using satellite data
    plt.plot(errors['NRSNSI_diff_mean'],np.arange(0,140),linewidth=4,
            label='NRSNSI (n prof = ' + str(errors['NRSNSI_n_profiles_used']) + ')')
    plt.plot(errors['NRSMAI_diff_mean'],np.arange(0,140),linewidth=4,
            label='NRSMAI (n prof = ' + str(errors['NRSMAI_n_profiles_used']) + ')')
    plt.plot(errors['NRSROT_diff_mean'],np.arange(0,140),linewidth=4,
            label='NRSROT (n prof = ' + str(errors['NRSROT_n_profiles_used']) + ')')
    plt.plot(errors['NRSKAI_diff_mean'],np.arange(0,140),linewidth=4,
            label='NRSKAI (n prof = ' + str(errors['NRSKAI_n_profiles_used']) + ')')
    plt.plot(errors['PH100_diff_mean'],np.arange(0,140),linewidth=4,
            label='PH100 (n prof = ' + str(errors['PH100_n_profiles_used']) + ')')
    
    # plot error assuming satellite data representative
    plt.plot(errors['NRSNSI_diff_mean'][0:22],np.arange(0,22),linewidth=2,linestyle='--', c='w')
    plt.plot(errors['NRSMAI_diff_mean'][0:21],np.arange(0,21),linewidth=2,linestyle='--', c='w')
    plt.plot(errors['NRSROT_diff_mean'][0:27],np.arange(0,27),linewidth=2,linestyle='--', c='w')
    plt.plot(errors['NRSKAI_diff_mean'][0:41],np.arange(0,41),linewidth=2,linestyle='--', c='w')
    plt.plot(errors['PH100_diff_mean'][0:15],np.arange(0,15),linewidth=2,linestyle='--', c='w')

    plt.scatter(np.interp(depths['NRSNSI'],np.arange(0,140),errors['NRSNSI_diff_mean']),
                depths['NRSNSI'],marker='s',c='k',s=30,zorder=10)
    plt.scatter(np.interp(depths['NRSNSI'],np.arange(0,140),errors['NRSNSI_diff_mean']),
                depths['NRSNSI'],marker='s',c='b',s=4,zorder=11)
    plt.scatter(np.interp(depths['NRSMAI'],np.arange(0,140),errors['NRSMAI_diff_mean']),
                depths['NRSMAI'],marker='s',c='k',s=30,zorder=10)
    plt.scatter(np.interp(depths['NRSMAI'],np.arange(0,140),errors['NRSMAI_diff_mean']),
                depths['NRSMAI'],marker='s',c='orange',s=4,zorder=11)
    plt.scatter(np.interp(depths['NRSROT'],np.arange(0,140),errors['NRSROT_diff_mean']),
                depths['NRSROT'],marker='s',c='k',s=30,zorder=10)
    plt.scatter(np.interp(depths['NRSROT'],np.arange(0,140),errors['NRSROT_diff_mean']),
                depths['NRSROT'],marker='s',c='green',s=4,zorder=11)
    plt.scatter(np.interp(depths['NRSKAI'],np.arange(0,140),errors['NRSKAI_diff_mean']),
                depths['NRSKAI'],marker='s',c='k',s=30,zorder=10)
    plt.scatter(np.interp(depths['NRSKAI'],np.arange(0,140),errors['NRSKAI_diff_mean']),
                depths['NRSKAI'],marker='s',c='r',s=4,zorder=11)
    plt.scatter(np.interp(depths['PH100'],np.arange(0,140),errors['PH100_diff_mean']),
                depths['PH100'],marker='s',c='k',s=30,zorder=10)
    plt.scatter(np.interp(depths['PH100'],np.arange(0,140),errors['PH100_diff_mean']),
                depths['PH100'],marker='s',c='purple',s=4,zorder=11)

    plt.legend(fontsize=12,frameon=False)
    plt.xlim([0,1.5])
    plt.ylim([-2,115])
    plt.gca().invert_yaxis()
    plt.grid()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Mean absolute temp. difference (mooring-CTD) [$^\circ$C]',fontsize=16)
    plt.ylabel('Depth [m]',fontsize=16)

    plt.tight_layout()

    # add legend items for std and satellite influenced portion
    plt.fill_betweenx([22, 24], 0.95, 1.03, 
                    color='k', alpha=0.1)
    plt.annotate(
        "Standard Dev.", 
        xy=(1.05,23), 
        xytext=(1.05, 24),
        fontsize=12,
        color='black'
    )
    plt.plot(np.linspace(0.95,1.03,10),np.ones(10)*27,linewidth=4, c='b')
    plt.plot(np.linspace(0.95,1.03,10),np.ones(10)*27,linewidth=2,linestyle='--', c='w')
    plt.annotate(
        "Extrap. /satellite", 
        xy=(1.05,28), 
        xytext=(1.05, 28),
        fontsize=12,
        color='black'
    )
    
    if season == 'all':
        plt.savefig(repo_path + '/Plots/' + 
                    'ProfileComparisons/Mooring_mean_diff_comparison.png',dpi=300)
    else:
        plt.savefig(repo_path + '/Plots/' + 
                            'ProfileComparisons/Mooring_mean_diff_comparison_' + season + '.png',dpi=300)
                
    plt.close()
    
plot_error(errors,'all')
plot_error(errors_summer,'summer')
plot_error(errors_winter,'winter')

# %% Create figure for mean error (per site plots)

colors = ['blue','orange','green','red','purple']
ls = ['NRSNSI','NRSMAI','NRSROT','NRSKAI','PH100']
satdepths = [22,21,27,41,15]

for n in range(len(ls)):

    s = ls[n]
    cl = colors[n]

    plt.figure(figsize=(8, 8))

    # Plot the shaded area representing ± standard deviation
    plt.fill_betweenx(np.arange(0,140), errors[s + '_diff_mean'] - errors[s + '_diff_std'], 
                    errors[s + '_diff_mean'] + errors[s + '_diff_std'], 
                    color=cl, alpha=0.1)

    # plot error using satellite data
    plt.plot(errors[s + '_diff_mean'],np.arange(0,140),linewidth=4, c=cl,   
            label=s + ' (n prof = ' + str(errors[s + '_n_profiles_used']) + ')')
  
    # plot error assuming satellite data representative
    plt.plot(errors[s + '_diff_mean'][0:satdepths[n]],np.arange(0,satdepths[n]),linewidth=2,linestyle='--', c='w')


    plt.scatter(np.interp(depths[s],np.arange(0,140),errors[s + '_diff_mean']),
                depths[s],marker='s',c='k',s=30,zorder=10)
    plt.scatter(np.interp(depths[s],np.arange(0,140),errors[s + '_diff_mean']),
                depths[s],marker='s',c=cl,s=4,zorder=11)


    plt.legend(fontsize=12,frameon=False)
    plt.xlim([0,1.5])
    plt.ylim([-2,115])
    plt.gca().invert_yaxis()
    plt.grid()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Mean absolute temp. difference (mooring-CTD) [$^\circ$C]',fontsize=16)
    plt.ylabel('Depth [m]',fontsize=16)
    plt.tight_layout()

    # add legend items for std and satellite influenced portion
    plt.fill_betweenx([4, 5], 0.97, 1.05, 
                    color=cl, alpha=0.1)
    plt.annotate(
        "Standard Dev.", 
        xy=(1.07,5.5), 
        xytext=(1.07, 5.5),
        fontsize=12,
        color='black'
    )
    plt.plot(np.linspace(0.97, 1.05,10),np.ones(10)*7,linewidth=4, c=cl)
    plt.plot(np.linspace(0.97, 1.05,10),np.ones(10)*7,linewidth=2,linestyle='--', c='w')
    plt.annotate(
        "Extrap. /satellite", 
        xy=(1.07,8), 
        xytext=(1.07, 8),
        fontsize=12,
        color='black'
    )

    plt.savefig(repo_path + '/Plots/ProfileComparisons/Mooring_mean_diff_comparison_' + s +'.png',dpi=300)
    plt.close()



# %% Compare satellite and CTD data
# sat SST vs near-surface CTD temp

for s in sites:
     
    print(s)
    
    # get mooring data
    Tm = pairs_ds[s + '_TEMP_mooring']
    Dm = pairs_ds[s + '_DEPTH_mooring']
    pnm = pairs_ds[s + '_profile_n_mooring']
    tm = pairs_ds[s + '_TIME_mooring']
    un_profs = np.int32(np.unique(pnm))
    # get CTD data
    CTDds = xr.open_dataset(repo_path + '/Data/' + 
                            s + '_CTD-M_pair.nc')
    TCTD = CTDds['TEMP_CTD'].values
    DCTD = CTDds['DEPTH_CTD'].values
    pnCTD = CTDds['profile_n_CTD'].values
    tCTD = CTDds['TIME_CTD'].values

    CTD_points = []
    CTD_depth = []
    sat_points = []

    for n in range(len(un_profs)):
        up = un_profs[n]
        print(up)
        # select data for profile
        cCTD = pnCTD == up;
        cM = pnm == up;
        c = np.argsort(Dm[cM])
        # select mooring and CTD data
        MD = Dm[cM][c]; MT = Tm[cM][c]
        CTDD = DCTD[cCTD]; CTDT = TCTD[cCTD]
        
        if CTDD.min() < 10:
            if np.isfinite(CTDT[CTDD == CTDD.min()][0]):
                CTD_points.append(CTDT[CTDD == CTDD.min()][0])
                CTD_depth.append(CTDD.min())
            else:
                CTD_points.append(CTDT[CTDD == CTDD.min()+1][0])
                CTD_depth.append(CTDD.min()+1)
            try:
                sat_points.append(MT[MD == 0.2][0])
            except:
                sat_points.append(np.nan)
                
    # create plot   
    plt.figure() 
    plt.scatter(CTD_points,sat_points,c=CTD_depth)
    plt.plot(np.arange(15,30),np.arange(15,30),c='black')  
    plt.colorbar(label='CTD depth (m)')
    plt.title(s + ' satellite vs CTD') 
    plt.grid('on') 
    plt.xlabel('CTD temperature [$^\circ$C]')
    plt.ylabel('Satellite SST [$^\circ$C]')
    
    # Calculate R2
    c = np.logical_and(np.isfinite(sat_points),
                       np.isfinite(CTD_points))
    r2 = r2_score(np.array(sat_points)[c], np.array(sat_points)[c])
    # Calculate RMSE
    rmse = np.sqrt(mean_squared_error(np.array(CTD_points)[c], np.array(sat_points)[c]))
        
    plt.savefig((repo_path + '/Plots/sat_vs_CTD/' + s + '.png'),dpi=300)
    plt.close('all')

# %% Plot mean CTD profiles for seasons

def calc_mean_profile(errors, sites):

    for s in sites:
        
        print(s)
        
        # get mooring data
        Tm = errors[s + '_mooring_temp']
        Dm = errors[s + '_mooring_depth']
        # get CTD data
        TCTD = errors[s + '_CTD_temp']
        DCTD = errors[s + '_CTD_depth']
        
        # get average binned mooring and CTD profiles
        D = []
        m_mean = []; m_std = [];
        CTD_mean = []; CTD_std = [];
        for d in range(0,141):
            c = np.logical_and(Dm > d,
                               Dm <= d+1)
            D.append(d+0.5)
            m_mean.append(np.nanmean(Tm[c]))
            m_std.append(np.nanstd(Tm[c]))
            CTD_mean.append(np.nanmean(TCTD[c]))
            CTD_std.append(np.nanstd(TCTD[c]))
        
        ds = {}
        ds[s + '_mooring_mean'] = np.array(m_mean)
        ds[s + '_mooring_std'] = np.array(m_std)
        ds[s + '_CTD_mean'] = np.array(CTD_mean)
        ds[s + '_CTD_std'] = np.array(CTD_std)
        ds[s + '_depth'] = np.array(D)
        
    return ds

means = calc_mean_profile(errors,sites)
means_summer = calc_mean_profile(errors_summer,sites)
means_winter = calc_mean_profile(errors_winter,sites)

# %%
