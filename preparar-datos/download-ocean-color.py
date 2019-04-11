import pandas as pd
from netCDF4 import Dataset
import datetime
import os
import math
import numpy as np
import requests
import json 

errors = []
meta = {}
lon_range = [-123,-87] # [-92.32910156250001, -90.6976318359375]
lat_range = [12,37.5] # [ 13.159725022841753,  14.689881366618774]
# prefix = "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/V"
# suffix = ".L3m_DAY_JPSS1_CHL_chlor_a_4km.nc"
def fetchData(file_prefix, file_suffix, query_date):
    yday = query_date.timetuple().tm_yday
    yearday = "{year}{day:03d}".format(year=query_date.year, day=yday)
    url = file_prefix + yearday + file_suffix
    print("Downloading " + url)
    os.system(" ".join(['wget', url,
                                "-O", "/rawdata/data.nc"
            ]) )
    if os.path.exists("/rawdata/data.nc"):
        try:
            dataset = Dataset("/rawdata/data.nc", "r", format="NETCDF4")
        except: 
            print("ERROR: Invalid data.")
            errors.append({
                "prefix": file_prefix,
                "suffix": file_suffix,
                "date": query_date
            })
            return None
    else:
        print("ERROR: Data file could not be downloaded")
        errors.append({
                "prefix": file_prefix,
                "suffix": file_suffix,
                "date": query_date
            })
        return None
    lons = np.argwhere((dataset.variables["lon"][:]>lon_range[0]) & (dataset.variables["lon"][:]<lon_range[1]))
    lats = np.argwhere((dataset.variables["lat"][:]>lat_range[0]) & (dataset.variables["lat"][:]<lat_range[1]))
    key = list(dataset.variables.keys())[0]
    subset = dataset.variables[key][min(lats)[0]:max(lats)[0],min(lons)[0]:max(lons)[0]]
    
    print("Data shape: ", subset.shape, " data var: ", key)
    alons, alats = dataset.variables["lon"][min(lons)[0]:max(lons)[0]], dataset.variables["lat"][min(lats)[0]:max(lats)[0]]

    lons, lats = np.meshgrid(alons, alats)
    adataframe = pd.DataFrame(data = { "value": subset.flatten(), "lons": lons.flatten(), "lats": lats.flatten() })
    yearday = query_date.year*1000 + yday
    adataframe["YearDay"] =  yearday
    instr_abbrev = dataset.instrument[0:3] + dataset.platform[0]
    adataframe["Instrument"] = instr_abbrev
    if instr_abbrev not in meta:
        meta[instr_abbrev] = {}
    if key not in meta[instr_abbrev]:
        meta[instr_abbrev][key] = {}
    meta[instr_abbrev][key][yearday] = {
            "latitude_step" : dataset.latitude_step,
            "longitude_step": dataset.longitude_step
        }
    dataset.close()
    os.system(" ".join(['rm', "/rawdata/data.nc"]) )
    
    return adataframe

# MODIS Aqua and Terra:
data_catalog = {
        "MODA": [
                {
                    "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A",
                    "suffix": ".L3m_DAY_CHL_chlor_a_4km.nc"
                },
                {
                    "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A",
                    "suffix": ".L3m_DAY_SST_sst_4km.nc"
                },
                {
                    "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A",
                    "suffix": ".L3m_DAY_GSM_bbp_443_gsm_4km.nc"
                }
            ] 
    }

# Lets download only 2018
d = datetime.date(2018,1,1)
df = datetime.date(2018,1,5)
dd = datetime.timedelta(days=1)
data = pd.DataFrame()
while d <= df:
    for mission, variables in data_catalog.items():
        for variable in variables:
            subdata = fetchData(variable["prefix"], variable["suffix"], d)
            data = pd.concat([data, subdata], ignore_index=True)
    d = d + dd
