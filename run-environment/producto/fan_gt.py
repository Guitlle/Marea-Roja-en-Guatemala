import math
import numpy as np
import requests
from PIL import Image
from io import BytesIO
import matplotlib as mlp
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import datetime
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



def getImageCluster(lat_deg, lon_deg, delta_lat,  delta_long, zoom):
    smurl = r"http://a.tile.openstreetmap.org/{0}/{1}/{2}.png"
    xmin, ymax =deg2num(lat_deg, lon_deg, zoom)
    xmax, ymin =deg2num(lat_deg + delta_lat, lon_deg + delta_long, zoom)

    Cluster = Image.new('RGB',((xmax-xmin+1)*256-1,(ymax-ymin+1)*256-1) ) 
    for xtile in range(xmin, xmax+1):
        for ytile in range(ymin,  ymax+1):
            imgurl=smurl.format(zoom, xtile, ytile)
            # print("Opening: " + imgurl)
            imgstr = requests.get(imgurl).content
            tile = Image.open(BytesIO(imgstr))
            Cluster.paste(tile, box=((xtile-xmin)*256 ,  (ytile-ymin)*255))

    return Cluster, [num2deg(xmin, ymax+1, zoom), num2deg(xmax+1, ymin, zoom)]

str_date = os.environ.get("DATE")

if str_date is None:
    query_date = datetime.datetime.today()
else:
    query_date = datetime.datetime.strptime(str_date, "%Y-%m-%d")                                                                                               

yday = query_date.timetuple().tm_yday
yearday = "{year}{day:03d}".format(year=query_date.year, day=yday)

os.system(" ".join(['wget', "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/V" + yearday + ".L3m_DAY_SNPP_CHL_chlor_a_4km.nc",
                               "-O", "/rawdata/data.nc"
           ]) )

d = Dataset("/rawdata/data.nc", "r", format="NETCDF4")


import geopandas as gp

lon_range = [-94,-87] # [-92.32910156250001, -90.6976318359375]
lat_range = [12,18] # [ 13.159725022841753,  14.689881366618774]
lons = np.argwhere((d.variables["lon"][:]>lon_range[0]) & (d.variables["lon"][:]<lon_range[1]))
lats = np.argwhere((d.variables["lat"][:]>lat_range[0]) & (d.variables["lat"][:]<lat_range[1]))
gtsubset = d.variables["chlor_a"][min(lats)[0]:max(lats)[0],min(lons)[0]:max(lons)[0]]
print(gtsubset.shape)

mapbg, bbox = getImageCluster(lat_range[0], lon_range[0], lat_range[1]-lat_range[0], lon_range[1]-lon_range[0], 7)

mlp.rcParams["figure.figsize"] = (18, 18)
mlp.rcParams["grid.alpha"] = 0.4
fig,ax = plt.subplots()
plt.xlim(lon_range[0], lon_range[1])
plt.ylim(lat_range[0], lat_range[1])
cmap = mlp.cm.jet
cmap.set_bad('black',0.3)
ax.imshow(mapbg, extent=(bbox[0][1], bbox[1][1], bbox[0][0], bbox[1][0]), interpolation = "gaussian")
satimg = ax.imshow(gtsubset, cmap = cmap, extent=(lon_range[0], lon_range[1], lat_range[0], lat_range[1]), 
                    norm=mlp.colors.LogNorm(vmin=gtsubset.min(), vmax=gtsubset.max()) )
plt.colorbar(satimg, shrink=0.5)

plt.savefig("/output/SNPP_chlor_a_gtm_"+yearday+".png")
