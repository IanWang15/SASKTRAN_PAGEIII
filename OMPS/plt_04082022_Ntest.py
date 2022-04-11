import sasktran as sk
import numpy as np
import matplotlib.pyplot as plt


tan_alts_km = np.arange(0.5, 60, 0.25)

geo = sk.VerticalImage()
geo.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tan_alts_km, mjd=54372, locallook=0)

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

# Grid to define our 2D cloud on
atmosphere_alts_m = np.arange(0, 100001, 1000)

# Angles in degrees relative to the TP, 0=TP, negative angles towards the observer
atmosphere_angles_deg = np.arange(-5, 5.01, 0.5)

cloud_numden = np.zeros((len(atmosphere_angles_deg), len(atmosphere_alts_m)))
print('cloud_numden ',np.shape(cloud_numden))

nclist = [0,0.005,0.01,0.05,0.1]

fig, ax = plt.subplots(figsize=(3, 7))
for inc in nclist:
# Localize the cloud +/1 degree from the TP from 10-12 km
    cloud_numden[9:11, 10:12] = inc

# Create the cloud climatology
    cloud_climatology = sk.ClimatologyUserDefined2D(atmosphere_angles_deg, atmosphere_alts_m,
                                                {'cloud': cloud_numden},
                                                normal_vector=normal,
                                                reference_vector=reference)

    atmo['cloud'] = sk.Species(cloud_opt_prop, cloud_climatology)

# Create the engine and do the calculation
    engine = sk.EngineHR(geometry=geo, atmosphere=atmo, wavelengths=[750])

# 2D atmosphere
    engine.atmosphere_dimensions = 2

    engine.configure_for_cloud(500, atmosphere_alts_m)

# Match the engine grid with our cloud grid
    engine.options['opticalnormalandreference'] = normalandreference
    engine.options['opticalanglegrid'] = atmosphere_angles_deg

    radiance = engine.calculate_radiance()

    if inc == 0:
        ax.plot(radiance[0], tan_alts_km, label='cloud free')
    else:
        ax.plot(radiance[0], tan_alts_km, label='N: '+str(inc))

#    if inc == 0:
#        ax.plot(cloud_numden[10,:], atmosphere_alts_m/1000., '--', label='cloud free')
#    else:
#        ax.plot(cloud_numden[10,:], atmosphere_alts_m/1000., '--', label='N: '+str(inc))

ax.set_xlabel('radiance')
ax.set_ylabel('altitude (km)')
ax.set_ylim(0,20.1)
    #ax.set_xlim(-0.005,0.121)
#ax.set_xticks([0, 0.04, 0.08, 0.12])
    #ax.set_title('N: '+str(inc)+'cm^-3')
plt.legend(loc='upper right')

pngname = '../fig/'+'fig_cloud_2d_N.png'
plt.savefig(pngname, bbox_inches='tight')

