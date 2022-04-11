import numpy as np
from pyhdf.SD import SD,SDC
import matplotlib.pyplot as plt

def statsnp(var):
    a = np.nanmax(var)
    b = np.nanmin(var)
    c = np.nanmedian(var)
    d = np.nanstd(var)
    e = np.nanmean(var)
    f = (c**2. + d**2.)
    print('max: ',a,'min: ', b,'median: ', c,'std: ', d, 'mean: ', e, 'uncertainty: ', f)

file00 = 'CAL/CAL_LID_L1-Standard-V4-10.2012-08-16T18-51-03ZD.hdf'

dfdir = '../dat/'+file00
#hf00 = SD(dfdir,SDC.WRITE)
hf00 = SD(dfdir,SDC.READ)
datasets_dic = hf00.datasets()
#for idx,sds in enumerate(datasets_dic.keys()):
#    print(idx,sds)

def readvar(varname):
    i00_obj = hf00.select(varname)
    i00 = i00_obj.get()

    add_offset,scale_factor = 0., 0.
    for key, value in i00_obj.attributes().items():
        if key == 'radiance_offsets':#,'add_offset':
            add_offset = value
        if key == 'radiance_scales':#,'scale_factor':
            scale_factor = value

#    varfillvalue = i00_obj.getfillvalue()
    i00_obj.endaccess()

    return i00#,add_offset,scale_factor,varfillvalue

# read variable
#var532 = readvar('Integrated_Attenuated_Backscatter_532')
var1064 = readvar('Attenuated_Backscatter_1064')
#vartopalt = readvar('Layer_Top_Altitude')
#varbasealt = readvar('Layer_Base_Altitude')
lat = readvar('Latitude')
lon = readvar('Longitude')

statsnp(lat)

statsnp(lon)

print(np.shape(var1064))
#print(np.shape(vartopalt))
print(np.shape(lat))
#print(np.shape(lon))


import h5py
import numpy as np
import matplotlib.pyplot as plt

# the data from OMPS
#file00 = 'OMPS-NPP_LP-L1G-EV_v2.5_2014m0219t072848_o11988_2016m0626t235315.h5'
file00 = 'OMPS/OMPS-NPP_LP-L1G-EV_v2.5_2012m0816t183034_o04163_2016m0622t105318.h5'

filename = '../dat/'+file00

#ds_disk = xr.open_dataset(filename)
#print(ds_disk)

with h5py.File(filename, 'r') as data:
#    for group in data.keys() :
#        print (group)

    ds_date = data['GRIDDED_DATA']['Date'][()]
    ds_time = data['GRIDDED_DATA']['Time'][()]
    ds_lat = data['GEOLOCATION_DATA']['Latitude_35km'][()]
    ds_lon = data['GEOLOCATION_DATA']['Longitude_35km'][()]


ntime = 179 # OMPS event index
slit = 1 # central slit

print(ds_lat[ntime,slit])
print(ds_lon[ntime,slit])

ompslat = ds_lat[:,slit]
ompslon = ds_lon[:,slit]


# ---------------------------------------------------------------------
import cartopy.crs as ccrs
#ds_disk = xr.open_dataset("saved_on_disk.nc")
import glob

data = np.ones(len(lat))
ompsdata = np.ones(len(ompslat)) + 1

fig = plt.figure(figsize=(6, 6))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-101, -69, 4, 31], crs=ccrs.PlateCarree())
#             west , east, south, north

mm = ax.scatter(lon, lat, c=data, alpha = 0.5, cmap='rainbow',s=15)
ax.scatter(ompslon, ompslat, c=ompsdata, alpha = 0.5, cmap='rainbow',s=15)
ax.coastlines()
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

cbar_ax = fig.add_axes([0.98, 0.25, 0.03, 0.5])
#cbar_ax.set_label('Cloud Top Temperature (K)')
plt.colorbar(mm, cax=cbar_ax)

#ax.set_global()

print('start saving')

pngname = "../fig/"+"fig_glob_orbit.png"
print("save ", pngname)
plt.savefig(pngname, bbox_inches='tight')
plt.show()
