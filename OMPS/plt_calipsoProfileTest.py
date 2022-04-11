import numpy as np
from pyhdf.SD import SD,SDC
from matplotlib.colors import LogNorm

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
var532 = readvar('Total_Attenuated_Backscatter_532')
var1064 = readvar('Attenuated_Backscatter_1064')
#vartopalt = readvar('Layer_Top_Altitude')
#varbasealt = readvar('Layer_Base_Altitude')
lat = readvar('Latitude')
lon = readvar('Longitude')
time = readvar('Profile_UTC_Time')

statsnp(var532)

statsnp(var1064)

print(np.shape(var1064))
#print(np.shape(vartopalt))
print(np.shape(lat))
#print(np.shape(lon))

#import pdb
#pdb.set_trace() 
#print('stop here')


latest = lat[lat[:,0].searchsorted(-52.)]
idxlat = np.where(lat[:,0] == latest)
print(idxlat)

latest = lat[lat[:,0].searchsorted(-71.)]
idxlat = np.where(lat[:,0] == latest)
print(idxlat)

latest = lat[lat[:,0].searchsorted(-31.)]
idxlat = np.where(lat[:,0] == latest)
print(idxlat)

#lonest = lon[lon[:,0].searchsorted(83.364)]
#idxlon = np.where(lon[:,0] == lonest)
print(idxlat)
#print(idxlon)
print(lat[idxlat,0])
#print(lon[idxlon,0])

loc0 = 31000

loc1 = int(np.squeeze(idxlat)-30)
loc2 = int(np.squeeze(idxlat)-20)
print('lat ',lat[loc0,0])
print('lon ',lon[loc0,0])
print('time',time[loc0,0])

#print('11250 ',np.argmax(var532[11250,:]))
#print('11750 ',np.argmax(var532[11750,:]))
#print('11000 ',np.argmax(var532[11000,:]))
#print('12010 ',np.argmax(var532[12010,:]))
#print('time ',time[11000],time[12010])

import matplotlib.pyplot as plt

imgmask = np.ma.array(var532,mask=var532 < -1)
img = np.rot90(imgmask[loc0:loc0+1010,50:])
img = np.rot90(img)
img = np.rot90(img)
#img.clip(min=0)

#img0=np.zeros(shape=(np.shape(img)[0],np.shape(img)[1]))

fig, ax = plt.subplots(figsize=(14, 8))
im = ax.imshow(img, cmap='rainbow', interpolation='nearest', norm=LogNorm(vmin=0.0001, vmax=0.1))
# cmap=plt.cm.viridis
ax.set_xticks([1010-0, 1010-500, 1010-1000])
ax.set_xticklabels([str(lat[loc0,0].round(2))+', '+str(lon[loc0,0].round(2)),str(lat[loc0+500,0].round(2))+', '+str(lon[loc0+500,0].round(2)),\
    str(lat[loc0+1000,0].round(2))+', '+str(lon[loc0+1000,0].round(2))], fontsize=18)
ax.set_yticks([533-0, 533-100, 533-200, 533-300, 533-400, 533-500])
ax.set_yticklabels(['0','6','12','18','24','30'], fontsize=18)
ax.set_xlabel('latitude (degree), longitude (degree)', fontsize=18)
ax.set_ylabel('altitude (km)', fontsize=18)

cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])
cbar = plt.colorbar(im, cax=cbar_ax)

pngname = "../fig/"+"fig_imshow_simple"+str(loc0)+".png"
plt.savefig(pngname, bbox_inches='tight')
print(pngname+' is saved')
