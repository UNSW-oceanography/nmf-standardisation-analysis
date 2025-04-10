# Vertical Mooring Configuration Analysis
## IMOS National Mooring Facility

This repository contains code to identify depths in the water column where sensors could potentially be deployed in future configurations. The code compares CTD profiles with interpolated profiles using mooring measurements at discrete depths, and estimates the difference between them. 

The code can be applied to data collected at any of the NMF mooring sites, but the methodology is heavily-dependent on the number of CTD profiles available at a site. 

## Repository Structure

### Code/

<table style="border: 4px solid black; border-collapse: collapse;">
  <thead>
    <tr>
      <th style="border: 4px solid black; padding: 8px;">Script Name</th>
      <th style="border: 4px solid black; padding: 8px;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">CTD-Mooring _Comparisons.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        Use this script to explore data availability at the sites, estimate the mixed layer depth, get the mean error between interpolated mooring profiles and CTD profiles, and create plots.  
      </td>
    </tr>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">CTD-mooring-pairing.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        This script matches the mooring data with the CTD profile close in time and saves the matched data in a NetCDF file for each mooring site. It also does some additional QC and matches SST data close in time and space to the CTD-mooring pair. 
      </td>
    </tr>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">getTEMP_LTSPs.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        Loads the latest long time series product for each mooring site and save them as xarray data sets inside a zipped pickle called 'TEMP_LTSPs.pkl.gz'. This script can take some time to complete. 
      </td>
    </tr>
  </tbody>
</table>

### Code/util

<table style="border: 4px solid black; border-collapse: collapse;">
  <thead>
    <tr>
      <th style="border: 4px solid black; padding: 8px;">Script Name</th>
      <th style="border: 4px solid black; padding: 8px;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">PickleStuff.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        Functions to save/load pickle files.
      </td>
    </tr>
  </tbody>
</table>

### Code/util/CTD_AGG

Use this code to create an aggregated CTD profile product.

<table style="border: 4px solid black; border-collapse: collapse;">
  <thead>
    <tr>
      <th style="border: 4px solid black; padding: 8px;">Script Name</th>
      <th style="border: 4px solid black; padding: 8px;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">GetCTDdata.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        Run this script to create the aggregated CTD profile products for each site defined in the `sites.yaml` config file. A dictionary called 'CTDdata' is created to store the products for each site, and later saved as a pickle file. 
      </td>
    </tr>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">aggregated_profiles.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        The functions used to create the aggregated CTD profile product. Contains a function to gather the thredds URLs for each individual CTD profile for a specific mooring site, and another function to aggregate those profiles into a dictionary called 'AggProfs'. 
      </td>
    </tr>
    <tr>
      <td style="border: 4px solid black; padding: 8px;">PickleStuff.py</td>
      <td style="border: 4px solid black; padding: 8px;">
        Functions to save/load pickle files. 
      </td>
    </tr>
  </tbody>
</table>

### Code/util/holteandtalley

This folder contains the holteandtalley package that estimates the mixed layer depth for a profile using various methodologies. You can only pip install this package, and I was using conda at the time so I instead imported the package from this folder.


## Setup

Create a virtual environment and install the packages listed in the `requirements.txt` file. 

Update the `sites.yaml` to include the mooring sites that you want to analyse. Each key is the mooring site name, and there is a subkey for the coordinates and depths. The coordinates are available in the metadata of the long time series product NetCDF files. The depths were determined by looking at the histogram of all DEPTH collected and selecting peaks by eye. 

## Methodology (and scripts)

### 1. Create an aggregated CTD product
______________________________

Update the `sites.yaml` to include additional mooring sites if you need, or leave as is. 

Edit the `repo_path` variable in each script to the path of your local version of the repository. 

Run `GetCTDdata.py` to create aggregated CTD profile products at each of the mooring sites.  

### 2. Pair the CTD profiles with Mooring interpolated profiles
______________________________

Create a folder called 'Data' in the repository. 

Run `getTEMP_LTSPs.py` to create a pickle file containing the temperature long time series products for all of the mooring sites. 

Then run `CTD-mooring-pairing.py` to pair the aggregated CTD profile products with the LTSPs, add the satellite SST, and do some additional QC. 

### 3. Calculate 'error' and create plots
______________________________

Go through `CTD-Mooring_Comparisons.py` section-by-section. 

You will need to create a folder called `Plots` and sub-folders called `DataAvailability`, `ProfileComparisons` and `sat_vs_CTD`.

Some editing might be required in some later sections. 

## Next Steps (Following a discussion at the QAQC summit 2025)

- Seasonal analysis – how does the mean CTD-mooring difference vary by season?​ keep in mind the number of profiles available at a particular site for each season. 
- Difference between the MLD estimated using CTD vs using mooring interpolated profiles. Keep in mind that the estimated MLD for the profiles are not always correct.
- Final tweaking of the methodology, further validation, checks:
  - Further inspection of profiles
  - Fix binning method / QC at certain depths (investigate) - look at profile comparison plots
  - Explore sensitivity of the overall error profile to number of profile pairs (look at histogram of profile availability over depth, some bins need to be excluded?)
