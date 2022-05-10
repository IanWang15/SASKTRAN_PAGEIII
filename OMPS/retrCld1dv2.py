import h5py
import numpy as np
import matplotlib.pyplot as plt

file00 = 'OMPS-NPP_LP-L1G-EV_v2.5_2014m0219t072848_o11988_2016m0626t235315.h5'

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
    ds_alt = data['GRIDDED_DATA']['TangentHeight'][()]
# the following also works
#    varlat = data['GEOLOCATION_DATA/Latitude_35km'][...]

slit = 1 # central slit
wavelength1 = .675 # microns
wavelength2 = .868 # microns
ClosestWavelength1 = WavelengthGrid[WavelengthGrid.searchsorted(wavelength1)]
# Identify the wavelength in the superset closest to the target wavelength
wavelengthidx1 = np.where(WavelengthGrid == ClosestWavelength1)
# Locate the index in the WavelengthGrid array

ClosestWavelength2 = WavelengthGrid[WavelengthGrid.searchsorted(wavelength2)]
wavelengthidx2 = np.where(WavelengthGrid == ClosestWavelength2)

ntime = 100
radimg = ds_rad[ntime,slit,:,wavelengthidx2]
radimgMask = np.ma.array(radimg,mask=radimg < -998)

# The fill value used in this dataset is -999

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

initarr[radimgMask[0,0,:] < 0.01] = 0.

#print(initarr)

# ----------------------------------------------

import sasktran as sk
from sasktran.geometry import VerticalImage
from sasktran import OpticalProperty, Climatology

# altitude follows OMPS
altarr = np.arange(0.5,40.6,1.)
ncarr = initarr[:len(altarr)] * 0.5

altarrmeter = altarr * 1000.
#print (altarr)
#print(ncarr)

#print(altarrmeter)

radm = radimgMask[0,0,:][:len(altarr)]

# First recreate our geometry and atmosphere classes
geometry = VerticalImage()
geometry.from_sza_saa(sza=30, saa=60, lat=0, lon=0, tanalts_km=altarr, mjd=54372, locallook=0,
                      satalt_km=840)#, refalt_km=15)

atmosphere0 = sk.Atmosphere()

atmosphere0['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
atmosphere0.brdf = sk.Lambertian(0.3)
wavelengths = np.array([868.])

# And now make the engine
engine0 = sk.EngineHR(geometry=geometry, atmosphere=atmosphere0, wavelengths=wavelengths)

rad0 = engine0.calculate_radiance()

atmosphere = sk.Atmosphere()

# Turn off Rayleigh scattering and gaseous in the air
#atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
#atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
atmosphere.atmospheric_state = sk.MSIS90()
#atmosphere['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
#atmosphere['no2'] = sk.Species(sk.NO2Vandaele1998(), sk.Pratmo())

#configure_for_cloud(grid_spacing_m: float, cloud_altitudes: numpy.array, cloud_diffuse_spacing_m: float = 100, max_optical_depth_of_cell: float = 0.01, min_extinction_ratio_of_cell: float = 1)

size = 1.

particle_size_dist = sk.ClimatologyUserDefined(altitudes=[0,45000], \
              values={'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': [1.6, 1.6], \
              'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': [size, size]})

aerosol_optprop = sk.MieAerosol(particlesize_climatology=particle_size_dist, species='WATER')


# the aerosol layer between 14000 m and 15990 m
aerosol_clim = sk.ClimatologyUserDefined(altitudes=altarrmeter, \
                                    values={'SKCLIMATOLOGY_AEROSOL_CM3': ncarr})

atmosphere['aer'] = sk.Species(aerosol_optprop, aerosol_clim)


atmosphere.brdf = sk.Lambertian(0.3)

#wavenum = np.array([934., 935.])
#wavelengths = 1.0E7/wavenum
#wavelengths = np.array([868., 675.])
#wavelengths = np.array([868.])

# And now make the engine
engine1 = sk.EngineHR(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)

rad1 = engine1.calculate_radiance()

opt_prop = OpticalProperty('rayleigh')
#sigma = opt_prop.calculate_cross_sections(atmospheric_state=Climatology('msis90'), latitude=0, longitude=0, altitude=50000, mjd=54372,wavelengths = np.array([868.]))
atmospheric_state = Climatology('msis90')
sigma = opt_prop.calculate_cross_sections(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,wavelengths=[868.])
#geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)
print(sigma)

fig, ax = plt.subplots(figsize=(3, 7))
im = ax.plot(rad0[0,:],altarr)

ax.set_xlabel('radiance')
ax.set_ylabel('altitude (km)')
ax.set_title('cloud free')
ax.set_ylim(0,40)

pngname = "../fig/"+"fig_cloud_aerFree.png"
plt.savefig(pngname, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3, 7))
im = ax.plot(rad1[0,:],altarr)

ax.set_xlabel('radiance')
ax.set_ylabel('altitude (km)')
ax.set_title('iteration 0')
ax.set_ylim(0,40)

pngname = "../fig/"+"fig_cloud_iter0.png"
plt.savefig(pngname, bbox_inches='tight')

#print(rad0)
print(np.shape(rad0))
#print(radm)
print(np.shape(rad0[0,:]))
print(np.shape(radm))
#print(altarrmeter)
print(ncarr)

# ----------------------------------------
for itr in range(10):
    rad0 = rad0[0,:]
    im = radm[:]/radm[-1]
    ic = rad0[:]/rad0[-1]
    i1 = rad1[:]/rad1[-1]

    ym = (im - i1)/i1
    yc = (ic - i1)/i1

    ym = ym[0,:]
    yc = yc[0,:]

    print(np.shape(ym))
    print(np.shape(yc))

    for i in range(len(radm)):
        #print(abs(radm[i] - rad0[i]))
#        print(abs(ym[i] - yc[i]),yc[i])
        if (abs(ym[i] - yc[i])/yc[i] > 0.2) and (ym[i] > yc[i]) and (rad0[i] > 0):
            ncarr[i] = ncarr[i] * 1.5
        if (abs(ym[i] - yc[i])/yc[i] > 0.2) and (ym[i] < yc[i]) and (rad0[i] > 0):
            ncarr[i] = ncarr[i] / 2.

#    import pdb
#    pdb.set_trace()
#    print('stop here')
    print(ncarr)

    atmosphere0 = sk.Atmosphere()
    atmosphere0['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    atmosphere0.brdf = sk.Lambertian(0.3)
    wavelengths = np.array([868.])
    # And now make the engine
    engine0 = sk.EngineHR(geometry=geometry, atmosphere=atmosphere0, wavelengths=wavelengths)
    rad0 = engine0.calculate_radiance()

    atmosphere = sk.Atmosphere()
    atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    atmosphere.atmospheric_state = sk.MSIS90()

    size = 1.
    particle_size_dist = sk.ClimatologyUserDefined(altitudes=[0,45000], \
              values={'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': [1.6, 1.6], \
              'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': [size, size]})
    aerosol_optprop = sk.MieAerosol(particlesize_climatology=particle_size_dist, species='WATER')
    aerosol_clim = sk.ClimatologyUserDefined(altitudes=altarrmeter, \
                                    values={'SKCLIMATOLOGY_AEROSOL_CM3': ncarr})

    atmosphere['aer'] = sk.Species(aerosol_optprop, aerosol_clim)
    atmosphere.brdf = sk.Lambertian(0.3)
    wavelengths = np.array([868.])
    # And now make the engine
    engine1 = sk.EngineHR(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)
    rad1 = engine1.calculate_radiance()

#    print(rad0)


    fig, ax = plt.subplots(figsize=(3, 7))
    im = ax.plot(rad1[0,:],altarr)

    ax.set_xlabel('radiance')
    ax.set_ylabel('altitude (km)')
    ax.set_title('iteration ' + str(itr+1))
    ax.set_ylim(0,40)

    pngname = "../fig/"+"fig_cloud_iter"+str(itr+1)+".png"
    plt.savefig(pngname, bbox_inches='tight')

print(radm)

fig, ax = plt.subplots(figsize=(3, 7))
im = ax.plot(radm,altarr)

ax.set_xlabel('radiance')
ax.set_ylabel('altitude (km)')
ax.set_title('measurement')
ax.set_ylim(0,40)

pngname = "../fig/"+"fig_cloud_measurement.png"
plt.savefig(pngname, bbox_inches='tight')
