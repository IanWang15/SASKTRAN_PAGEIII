import sasktran as sk
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
    WavelengthGrid = data['/GRIDDED_DATA/WavelengthGrid'][...]

    ds_rad = data['GRIDDED_DATA']['Radiance'][()]
    ds_rflc = data['GRIDDED_DATA']['Reflectance'][()]
    ds_lat = data['GEOLOCATION_DATA']['Latitude_35km'][()]
    ds_lon = data['GEOLOCATION_DATA']['Longitude_35km'][()]
    ds_sza = data['GEOLOCATION_DATA']['SolarZenithAngle_35km'][()]
    ds_saa = data['GEOLOCATION_DATA']['SolarAzimuth_35km'][()]
    ds_salt = data['GEOLOCATION_DATA']['SpacecraftAltitude'][()]
    ds_alt = data['GRIDDED_DATA']['TangentHeight'][()]
    ds_date = data['GRIDDED_DATA']['Date'][()]
    ds_time = data['GRIDDED_DATA']['Time'][()]

# the following also works
#    varlat = data['GEOLOCATION_DATA/Latitude_35km'][...]

ntime = 87 # OMPS event index
slit = 1 # central slit
wavelength1 = .675 # microns
wavelength2 = .868 # microns
sza = float(ds_sza[ntime,slit])
saa = float(ds_saa[ntime,slit])
lat = float(ds_lat[ntime,slit])
lon = float(ds_lon[ntime,slit])
salt = float(ds_salt[ntime,slit])

print('lat ',lat)
print('lon ',lon)

date = ds_date[ntime,slit]
time = ds_time[ntime,slit]
#print(date)
year = int(date/10000)
month = int((date - year*10000)/100)
day = date - year*10000 - month*100
#print(year, month, day)
#print(time)
hour = int(time/3600)
minute = int((time - hour * 3600)/60)
second = time - hour * 3600 - minute * 60
#print(hour, minute, second)

# importing pandas as pd
import pandas as pd
# Create the Timestamp object
ts = pd.Timestamp(year = year,  month = month, day = day,
                  hour = hour, minute = minute, second = second, tz = 'UTC')
#print(ts) # Print the Timestamp object
# convert to julian date
jd = ts.to_julian_date()
# convert to modified julian date
mjd = jd - 2400000.5
#print(mjd)

ClosestWavelength1 = WavelengthGrid[WavelengthGrid.searchsorted(wavelength1)]
# Identify the wavelength in the superset closest to the target wavelength
wavelengthidx1 = np.where(WavelengthGrid == ClosestWavelength1)
# Locate the index in the WavelengthGrid array

ClosestWavelength2 = WavelengthGrid[WavelengthGrid.searchsorted(wavelength2)]
wavelengthidx2 = np.where(WavelengthGrid == ClosestWavelength2)

radimg = ds_rad[ntime,slit,:,wavelengthidx2]
# times, slit, tangent height, wavelength
# extract single event, central slit, wavelength 2: 0.868

# mask: negative values -> null
# The fill value used in this dataset is -999
radimgMask = np.ma.array(radimg,mask=radimg < -998)

print(ds_rad.shape)
# (180,3,101,300)
# times, slit, tangent height, wavelength
print(ds_alt.shape)
print(ds_lat.shape)

#print(ragimgMask)
#print(ds_alt[100,slit,:])

#print(radimgMask[0,0,:])
#print(np.shape(radimgMask[0,0,:]))

initarr = np.ones(len(radimgMask[0,0,:]))

initarr[radimgMask[0,0,:] < 0.03] = 0.

print('radimgMask',np.shape(radimgMask[0,0,:]))
print('initarr',np.shape(initarr))

# ----------------------------------------------
# The initial RT calculation with clear sky profile

#tan_alts_km = np.arange(0.5, 60, 0.25)
tan_alts_km = np.arange(0.5, 40.6, 1.)
print('tan_alts_km ',np.shape(tan_alts_km))

radm = radimgMask[0,0,:][:len(tan_alts_km)] # radiance from measurement
radm = np.squeeze(radm)

print(sza,saa,lat,lon,salt)
geo = sk.VerticalImage()
geo.from_sza_saa(sza=sza, saa=saa, lat=lat, lon=lon, tanalts_km=tan_alts_km, mjd=mjd, locallook=0, satalt_km=salt, refalt_km=35)

atmo = sk.Atmosphere()
atmo['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

# Create the cloud optical property
cloud_opt_prop = sk.BaumIceCrystal(50)

# Need two vectors to define the 2D plane for the cloud
# Unit vector at the tangent point (TP) is one of them
reference = geo.reftp / np.linalg.norm(geo.reftp)
print('geo.reftp ',geo.reftp)
print('reference',)

# The other is a normal vector to the line of sight plane
normal = np.cross(geo.reftp, geo.reflook)
normal /= np.linalg.norm(normal)

print('geo.reflook ',geo.reflook)
print('normal ', normal)

# SASKTRAN takes them in as a stacked vector
normalandreference = np.concatenate((normal, reference))

print('normalandreference',normalandreference)
#import pdb
#pdb.set_trace()
#print('stop here')

# Grid to define our 2D cloud on
atmosphere_alts_m = np.arange(0, 25001, 1000)

# Angles in degrees relative to the TP, 0=TP, negative angles towards the observer
atmosphere_angles_deg = np.arange(-5, 5.01, 0.5)

cloud_numden = np.zeros((len(atmosphere_angles_deg), len(atmosphere_alts_m)))
print('cloud_numden ',np.shape(cloud_numden))

cloud_climatology = sk.ClimatologyUserDefined2D(atmosphere_angles_deg, atmosphere_alts_m,
                                                {'cloud': cloud_numden},
                                                normal_vector=normal,
                                                reference_vector=reference)

atmo['cloud'] = sk.Species(cloud_opt_prop, cloud_climatology)

# Create the engine and do the calculation
engine = sk.EngineHR(geometry=geo, atmosphere=atmo, wavelengths=[868])

# 2D atmosphere
engine.atmosphere_dimensions = 2

engine.configure_for_cloud(500, atmosphere_alts_m)

# Match the engine grid with our cloud grid
engine.options['opticalnormalandreference'] = normalandreference
engine.options['opticalanglegrid'] = atmosphere_angles_deg

rad0 = engine.calculate_radiance() # radiance from clear sky profile

if (rad0.ndim == 2):
    rad0 = rad0[0,:]

# ----------------------------------------------------------------------------------
# first iteration

# number concentraion profile for RT calculation
ncarr = initarr[:len(np.squeeze(cloud_numden[1,:]))] * 0.01
#print(ncarr)



# Localize the cloud +/1 degree from the TP from 10-12 km
cloud_numden[9:11, :] = np.squeeze(ncarr)

# Create the cloud climatology
cloud_climatology = sk.ClimatologyUserDefined2D(atmosphere_angles_deg, atmosphere_alts_m,
                                                {'cloud': cloud_numden},
                                                normal_vector=normal,
                                                reference_vector=reference)

atmo['cloud'] = sk.Species(cloud_opt_prop, cloud_climatology)

# Create the engine and do the calculation
engine = sk.EngineHR(geometry=geo, atmosphere=atmo, wavelengths=[868])

# 2D atmosphere
engine.atmosphere_dimensions = 2

engine.configure_for_cloud(500, atmosphere_alts_m)

# Match the engine grid with our cloud grid
engine.options['opticalnormalandreference'] = normalandreference
engine.options['opticalanglegrid'] = atmosphere_angles_deg

rad1 = engine.calculate_radiance()
if (rad1.ndim == 2):
    rad1 = rad1[0,:]

#xs = cloud_opt_prop.calculate_cross_sections(sk.MSIS90(),latitude=lat, longitude=lon, altitude=tan_alts_km, mjd=mjd, wavelengths=868).total
#print(np.shape(xs))
#print(xs)

fig, ax1 = plt.subplots(figsize=(3, 7))
#ax[1].plot(rad1, tan_alts_km, label='1st', alpha=0.1)
#ax[0].plot(cloud_numden[10,:], atmosphere_alts_m/1000., label='iter: 1', alpha=0.1)

# --------------------------------------------------------------------------
# start iteration

for itr in range(6):
    print('iteratrion: ', itr)
    print('number concentration ', cloud_numden[10,:])

    im = radm[:]/radm[-1]
    ic = rad0[:]/rad0[-1]
    i1 = rad1[:]/rad1[-1]

    ym = (im - i1)/i1
    yc = (ic - i1)/i1

#    ym = ym[0,:]
#    yc = yc[0,:]

    print(np.shape(ym))
    print(np.shape(yc))

#    for i in range(len(ncarr)):
#        #print(abs(radm[i] - rad0[i]))
#        print(ym[i],yc[i],abs(ym[i] - yc[i])/yc[i])
#        if (abs(ym[i] - yc[i])/yc[i] > 1.2) and (ym[i] > yc[i]) and (rad0[i] > 0):
#            print(i,ym[i],yc[i],rad0[i],rad1[i],ncarr[i])
#            ncarr[i] = ncarr[i] * 3
#        if (abs(ym[i] - yc[i])/yc[i] > 1.2) and (ym[i] < yc[i]) and (rad0[i] > 0):
#            ncarr[i] = ncarr[i] / 5.

    for i in range(len(ncarr)):
        #print(abs(radm[i] - rad0[i]))
#        print(ym[i],yc[i],abs(ym[i] - yc[i])/yc[i])
        if (abs(radm[i] - rad1[i])/radm[i] > 0.1) and (radm[i] > rad1[i]) and (rad0[i] > 0):
#            print(i,ym[i],yc[i],rad0[i],rad1[i],ncarr[i])
            ncarr[i] = ncarr[i] * 1.2
#        if (abs(ym[i] - yc[i])/yc[i] > 1.2) and (ym[i] < yc[i]) and (rad0[i] > 0):
        if (abs(radm[i] - rad1[i])/radm[i] > 0.1) and (radm[i] < rad1[i]) and (rad0[i] > 0):
#            print(i,radm[i],rad0[i],rad1[i],ncarr[i])
            ncarr[i] = ncarr[i] / 1.3
        if (ncarr[i] > 0.3):
            ncarr[i] =0.3

# Localize the cloud +/1 degree from the TP from 10-12 km
    cloud_numden[9:11, :] = np.squeeze(ncarr)

# Create the cloud climatology
    cloud_climatology = sk.ClimatologyUserDefined2D(atmosphere_angles_deg, atmosphere_alts_m,
                                                {'cloud': cloud_numden},
                                                normal_vector=normal,
                                                reference_vector=reference)

    atmo['cloud'] = sk.Species(cloud_opt_prop, cloud_climatology)

# Create the engine and do the calculation
    engine = sk.EngineHR(geometry=geo, atmosphere=atmo, wavelengths=[868])

# 2D atmosphere
    engine.atmosphere_dimensions = 2

    engine.configure_for_cloud(500, atmosphere_alts_m)
# Match the engine grid with our cloud grid
    engine.options['opticalnormalandreference'] = normalandreference
    engine.options['opticalanglegrid'] = atmosphere_angles_deg

    rad1 = engine.calculate_radiance()

    if (rad1.ndim == 2):
        rad1 = rad1[0,:]

    print('radiance',np.shape(rad1))
    print(np.shape(tan_alts_km))

#plt.plot(radiance[0], tan_alts_km)

## repeat the same calculation with the cloud horizontally homogeneous
#cloud_homog = np.zeros_like(cloud_numden)
#cloud_homog[:, 10:12] = 0.01
#cloud_climatology['cloud'] = cloud_homog

#radiance = engine.calculate_radiance()

#    ax[1].plot(rad1, tan_alts_km, label='iter: '+str(itr+2), alpha=0.1+itr*0.05)
#    ax[0].plot(cloud_numden[10,:], atmosphere_alts_m/1000., label=str(itr+2), alpha=0.1+itr*0.05)

ax1.plot(rad1, tan_alts_km, color = 'red', label='iter time: '+str(itr+2), alpha=0.5)
ax2 = ax1.twiny()
ax2.plot(cloud_numden[10,:], atmosphere_alts_m/1000.,'--',color = 'k', label='N', alpha=0.5)

ax1.plot(rad0, tan_alts_km, color = 'gold', label='clear sky', alpha=0.5)
ax1.plot(radm, tan_alts_km, color = 'blue', label='observation', alpha=0.5)

ax2.set_xlabel('N (cm^-3)')
ax2.set_xlim(-0.002,0.031)
#ax[0,0].set_xticks([0, 0.004, 0.008, 0.012])

ax1.set_xlabel('radiance')
ax1.set_ylabel('altitude (km)')
ax1.set_ylim(0,30.1)
ax1.set_xlim(-0.005,0.165)
ax1.set_xticks([0, 0.04, 0.08, 0.12, 0.16])
    #ax.set_title('N: '+str(inc)+'cm^-3')
ax1.legend(loc='upper right', frameon=False)
ax2.legend(loc='upper left', frameon=False)

pngname = '../fig/'+'fig_cloud_2d_iter_'+str(ntime)+'.png'
plt.savefig(pngname, bbox_inches='tight')
