import math
import numpy as np
import requests
from PIL import Image
from io import BytesIO
import matplotlib as mlp
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import datetime
import pandas as pd
import os

# Helper functions to plot a OSM map as background
def deg2num(lat_deg, lon_deg, zoom):
  lat_rad = math.radians(lat_deg)
  n = 2.0 ** zoom
  xtile = int((lon_deg + 180.0) / 360.0 * n)
  ytile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
  return (xtile, ytile)

def num2deg(xtile, ytile, zoom):
  n = 2.0 ** zoom
  lon_deg = xtile / n * 360.0 - 180.0
  lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
  lat_deg = math.degrees(lat_rad)
  return (lat_deg, lon_deg)


def getImageClusterPreprocesado(lat_deg, lon_deg, delta_lat,  delta_long, zoom):
    # Esto carga una imagen preprocesada para el fondo del mapa que encaja con estas coordenadas:
    # lon_range = [-93.2,-87] 
    # lat_range = [12,18.5] 
    
    smurl = r"http://a.tile.openstreetmap.org/{0}/{1}/{2}.png"
    xmin, ymax =deg2num(lat_deg, lon_deg, zoom)
    xmax, ymin =deg2num(lat_deg + delta_lat, lon_deg + delta_long, zoom)
    
    Cluster = Image.open("app/fondo-gt-20190624.png")
    
    return Cluster, [num2deg(xmin, ymax+1, zoom), num2deg(xmax+1, ymin, zoom)]


def getImageCluster(lat_deg, lon_deg, delta_lat,  delta_long, zoom):
    smurl = r"http://a.tile.openstreetmap.org/{0}/{1}/{2}.png"
    xmin, ymax =deg2num(lat_deg, lon_deg, zoom)
    xmax, ymin =deg2num(lat_deg + delta_lat, lon_deg + delta_long, zoom)
    
    Cluster = Image.new('RGB',((xmax-xmin+1)*256-1,(ymax-ymin+1)*256-1) ) 
    for xtile in range(xmin, xmax+1):
        for ytile in range(ymin,  ymax+1):
            imgurl=smurl.format(zoom, xtile, ytile)
            # print("Opening: " + imgurl)
            imgstr = requests.get(imgurl, headers = {'User-Agent': 'Mozilla/5.0 (Windows; U; Windows NT 5.1; hu-HU; rv:1.7.8) Gecko/20050511 Firefox/1.0.4'}).content
            tile = Image.open(BytesIO(imgstr))
            Cluster.paste(tile, box=((xtile-xmin)*256 ,  (ytile-ymin)*255))

    return Cluster, [num2deg(xmin, ymax+1, zoom), num2deg(xmax+1, ymin, zoom)]


# Get date from environment variable DATE in the format 2019-01-01
str_date = os.environ.get("DATE")

if str_date is None:
    query_date = datetime.datetime.utcnow().date() - datetime.timedelta(days=1)
else:
    query_date = datetime.datetime.strptime(str_date, "%Y-%m-%d")                                                                                               

yday = query_date.timetuple().tm_yday
yearday = "{year}{day:03d}".format(year=query_date.year, day=yday)

import geopandas as gp

mlp.rcParams["font.family"] ='serif'
mlp.rcParams["figure.dpi"] = 175
mlp.rcParams["figure.figsize"] = (9, 9)
mlp.rcParams["grid.alpha"] = 0.5

def plot_gt(d):
    lon_range = [-93.2,-87] # [-92.32910156250001, -90.6976318359375]
    lat_range = [12,18.5] # [ 13.159725022841753,  14.689881366618774]
    lons = np.argwhere((d.variables["lon"][:]>lon_range[0]) & (d.variables["lon"][:]<lon_range[1]))
    lats = np.argwhere((d.variables["lat"][:]>lat_range[0]) & (d.variables["lat"][:]<lat_range[1]))
    gtsubset = d.variables["chlor_a"][min(lats)[0]:max(lats)[0],min(lons)[0]:max(lons)[0]]
    print(gtsubset.shape)

    mapbg, bbox = getImageClusterPreprocesado(lat_range[0], lon_range[0], lat_range[1]-lat_range[0], lon_range[1]-lon_range[0], 7)

    fig,ax = plt.subplots()
    plt.xlim(lon_range[0], lon_range[1])
    plt.ylim(lat_range[0], lat_range[1])
    cmap = mlp.cm.get_cmap("jet", 20)
    cmap.set_bad('black',0.3)
    norm = mlp.colors.LogNorm(vmin=0.01, vmax=40)  # (vmin=gtsubset.min(), vmax=gtsubset.max())
    ax.imshow(mapbg, extent=(bbox[0][1], bbox[1][1], bbox[0][0], bbox[1][0]), interpolation = "gaussian")
    satimg = ax.imshow(gtsubset, cmap = cmap, extent=(lon_range[0], lon_range[1], lat_range[0], lat_range[1]), 
                        norm=norm)
    cb = plt.colorbar(satimg, shrink=0.75)
    #cb.set_ticks([x+0.5 for x in norm.inverse([(x+0.0001) for x in range(0,30,2)]) ])
    #cb.set_ticklabels([x for x in range(0,30,2) ])
    cb.set_ticks([np.round(x,2) for x in norm.inverse([(x-0.5)/20 for x in range(1,21)])], update_ticks=True)
    cb.set_ticklabels([np.round(x,2) for x in norm.inverse([(x-0.5)/20 for x in range(1,21)])])
    cb.ax.tick_params(which = 'minor', length = 0)
    plt.grid(True)
    plt.tight_layout()
    return gtsubset, d.variables["lon"][min(lons)[0]:max(lons)[0]], d.variables["lat"][min(lats)[0]:max(lats)[0]]


# Fetch data from NASA's Oceancolor data repository
# VIIRS SNPP
os.system(" ".join(['wget', "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/V" + yearday + ".L3m_DAY_SNPP_CHL_chlor_a_4km.nc",
                               "-O", "/rawdata/data1.nc"
           ]) )

if os.path.exists("/rawdata/data1.nc"):
    try:
        d1 = Dataset("/rawdata/data1.nc", "r", format="NETCDF4")
    except: 
        print("ERROR: Invalid VIIRS SNPP data.")
    else:
        data, alons, alats = plot_gt(d1)
        lons, lats = np.meshgrid(alons, alats)
        plt.title("VIIRS-SNPP (chlor_a)\nClorofila (mg/m³) \n "+str(query_date.date()) )
        plt.savefig("/output/VIIRS-SNPP_chlor_a_gtm_"+str(query_date.date())+".png")
        pd.DataFrame(data = { "value": data.flatten(), "lons": lons.flatten(), "lats": lats.flatten() }).to_csv("/output/VIIRS-SNPP_chlor_a_gtm_data_"+str(query_date.date())+".csv")
        plt.clf()
        d1.close()
else:
    print("ERROR: VIIRS SNPP data could not be downloaded")

# VIIRS JPSS
os.system(" ".join(['wget', "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/V" + yearday + ".L3m_DAY_JPSS1_CHL_chlor_a_4km.nc",
                               "-O", "/rawdata/data2.nc"
           ]) )
if os.path.exists("/rawdata/data2.nc"):
    try:
        d2 = Dataset("/rawdata/data2.nc", "r", format="NETCDF4")
    except: 
        print("ERROR: Invalid VIIRS JPSS data.")
    else:
        data, alons, alats = plot_gt(d2)
        lons, lats = np.meshgrid(alons, alats)
        plt.title("VIIRS-JPSS (chlor_a)\nClorofila (mg/m³) \n "+str(query_date.date()) )
        plt.savefig("/output/VIIRS-JPSS1_chlor_a_gtm_"+str(query_date.date())+".png")
        pd.DataFrame(data = { "value": data.flatten(), "lons": lons.flatten(), "lats": lats.flatten() }).to_csv("/output/VIIRS-JPSS1_chlor_a_gtm_data_"+str(query_date.date())+".csv")
        plt.clf()
        d2.close()
else:
    print("ERROR: VIIRS JPSS data could not be downloaded")

# AQUA MODIS
os.system(" ".join(['wget', "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A" + yearday + ".L3m_DAY_CHL_chlor_a_4km.nc",
                               "-O", "/rawdata/data3.nc"
           ]) )

if os.path.exists("/rawdata/data3.nc"):
    try:
        d3 = Dataset("/rawdata/data3.nc", "r", format="NETCDF4")
    except: 
        print("ERROR: Invalid MODIS Aqua data.")
    else:
        data, alons, alats = plot_gt(d3)
        lons, lats = np.meshgrid(alons, alats)
        plt.title("MODIS Aqua (chlor_a)\nClorofila (mg/m³) \n "+str(query_date.date()) )
        plt.savefig("/output/MODIS-Aqua_chlor_a_gtm_"+str(query_date.date())+".png")
        pd.DataFrame(data = { "value": data.flatten(), "lons": lons.flatten(), "lats": lats.flatten() }).to_csv("/output/MODIS-Aqua_chlor_a_gtm_data_"+str(query_date.date())+".csv")
        plt.clf()
        d3.close()
else:
    print("ERROR: MODIS Aqua data could not be downloaded")

