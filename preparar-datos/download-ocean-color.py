import pandas as pd
from netCDF4 import Dataset
import datetime
import os
import sys
import math
import time
import numpy as np
import xarray
import requests
import json 

errors = []
meta = {}
# Incluir datos de california para estudiar los datos de SCOOS
#lon_range = [-123,-87] # [-92.32910156250001, -90.6976318359375]
#lat_range = [11.8,37.5] # [ 13.159725022841753,  14.689881366618774]
lon_range = [-93.2,-87] # [-92.32910156250001, -90.6976318359375]
lat_range = [12,18.5] # [ 13.159725022841753,  14.689881366618774]


rawdata_path = "./raw/"
output_path = "./output/"


# prefix = "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/V"
# suffix = ".L3m_DAY_JPSS1_CHL_chlor_a_4km.nc"
class OceanColorDataset:
    def __init__(self, file_prefix, file_suffix, query_date):
        self.file_prefix = file_prefix
        self.file_suffix = file_suffix
        self.query_date  = query_date
    
    def fetchData(self):
        yday = self.query_date.timetuple().tm_yday
        self.yearday = "{year}{day:03d}".format(year=self.query_date.year, day=yday)
        url = self.file_prefix + self.yearday + self.file_suffix
        print("Downloading " + url)
        os.system(" ".join(['wget', url,
                                    "-O", rawdata_path + "data.nc"
                ]) )
    def readData(self):
        global errors
        if os.path.exists(rawdata_path + "data.nc"):
            try:
                xrd = xarray.open_dataset(rawdata_path + "data.nc" , cache=False)
                print("QUERY DATE: ", self.query_date)
            except: 
                print("ERROR: Invalid data. ", self.query_date)
                self.dataset = None
                errors.append({
                    "prefix": self.file_prefix,
                    "suffix": self.file_suffix,
                    "date": self.query_date
                })
                
                return None
        else:
            print("ERROR: Data file could not be downloaded")
            self.dataset = None
            errors.append({
                    "prefix": self.file_prefix,
                    "suffix": self.file_suffix,
                    "date": self.query_date
                })
            return None
        self.dataset = xrd
    
    def save(self, file_prefix, lat_range, lon_range ):
        xrd = self.dataset.sel(lat = slice(lat_range[1],lat_range[0] ),lon=  slice(lon_range[0], lon_range[1]))

        xrd.to_netcdf(output_path + file_prefix + self.yearday + self.file_suffix)

    def close(self):
        self.dataset.close()
        os.system(" ".join(['rm', rawdata_path + "data.nc"]) )
        # Try reload the dataset just to force refresh (somehow netcdf4 does not read the actual file, it assumes it has the same data)
        try:
            xrd = xarray.open_dataset(rawdata_path + "data.nc" , cache=False)
            xrd.close()
        except: 
            pass
        
# MODIS Aqua and Terra:
data_catalog = {
        #"MODT": [
        #        {
        #            "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A",
        #            "suffix": ".L3m_DAY_CHL_chlor_a_4km.nc"
        #        },
        #        {
        #            "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A",
        #            "suffix": ".L3m_DAY_SST_sst_4km.nc"
        #        },
        #        {
        #            "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A",
        #            "suffix": ".L3m_DAY_GSM_bbp_443_gsm_4km.nc"
        #        }
        #    ],
        "MODT": [
                {
                    "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/T",
                    "suffix": ".L3m_DAY_CHL_chlor_a_4km.nc"
                },
                {
                    "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/T",
                    "suffix": ".L3m_DAY_SST_sst_4km.nc"
                },
                {
                    "prefix": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/T",
                    "suffix": ".L3m_DAY_GSM_bbp_443_gsm_4km.nc"
                }
            ]
    }


if __name__ == "__main__":
    if (len(sys.argv) == 3):
        d    = datetime.datetime.strptime(sys.argv[1], "%Y-%m-%d")  # datetime.date(2018,1,1)
        df   = datetime.datetime.strptime(sys.argv[2], "%Y-%m-%d")  # datetime.date(2018,1,5)
    else:
        print("Usage: \n  python download-ocean-color.py <date start> <date end>")
        exit()
    
    dd   = datetime.timedelta(days=1)
    data = pd.DataFrame()
    while d <= df:
        for mission, variables in data_catalog.items():
            for variable in variables:
                o = OceanColorDataset(variable["prefix"], variable["suffix"], d)                                                                                                                  
                o.fetchData()
                o.readData()
                if o.dataset is not None:
                    o.save(mission, lat_range, lon_range)
                time.sleep(1.5)
        d = d + dd
        time.sleep(2)
    print(errors)

# Make some quick tests
# Ignore this
if False:
    mission, variables =  data_catalog.items().__iter__().__next__()                                                                                                                  
    variable = variables[0]                                                                                                                                                           
    d = datetime.datetime(2018,1,15) 
    o = OceanColorDataset(variable["prefix"], variable["suffix"], d)                                                                                                                  
    o.fetchData()
    o.readData()
    o.save(mission, lat_range, lon_range)
    #o.dataset["chlor_a"].plot()
    print(o.dataset.attrs["time_coverage_end"])
    o.close()
    d = datetime.datetime(2018,12,13) 
    o = OceanColorDataset(variable["prefix"], variable["suffix"], d)                                                                                                                  
    o.fetchData()
    o.readData()
    o.save(mission, lat_range, lon_range)    #o.dataset["chlor_a"].plot()
    print(o.dataset.attrs["time_coverage_end"])
    o.close()
    
