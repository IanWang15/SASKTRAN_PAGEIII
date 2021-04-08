import sasktran as sk
from sasktran.geometry import VerticalImage
import numpy as np

# First recreate our geometry and atmosphere classes
geometry = VerticalImage()
geometry.from_sza_saa(sza=90, saa=60, lat=0, lon=0, tanalts_km=[5, 10, 15, 20], mjd=54372, locallook=0,
                      satalt_km=400, refalt_km=15)

atmosphere = sk.Atmosphere()

# Turn off Rayleigh scattering and gaseous in the air
#atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

#wavenum = np.array([934., 935.])
#wavelengths = 1.0E7/wavenum

#wavenum = np.array([550., 1022.])
wavelengths = np.array([1022., 550.])

# And now make the engine
engine0 = sk.EngineOCC(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)

tau0 = engine0.calculate_radiance()

# an aerosol layer
# height: between 14 km and 15 km
# effective radius = 0.1 um
# aerosol type: H2SO4

#f0 = open('./f_h2so4.txt','w')
earth = 6375.
l = np.sqrt((6375.+16.)**2 - (6375.+15.)**2)*2.

for i in range(50):
    size = np.round(0.05+0.05*i, 2)

    particle_size_dist = sk.ClimatologyUserDefined(altitudes=[0, 16000], \
                  values={'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': [1.6, 1.6], \
                  'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': [size, size]})

    aerosol_optprop = sk.MieAerosol(particlesize_climatology=particle_size_dist, species='H2SO4')

    # the aerosol layer between 10 m and 15990 m
    aerosol_clim = sk.ClimatologyUserDefined(altitudes=[0,10,15990,16000], \
                                        values={'SKCLIMATOLOGY_AEROSOL_CM3': [0,100,100,0]})

    atmosphere['aer'] = sk.Species(aerosol_optprop, aerosol_clim)

    # And now make the engine
    engine = sk.EngineOCC(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)

    tau = engine.calculate_radiance()

    #     size, color ratio (550 nm/1022 nm), extinction coefficient @ 1022 nm
    print(size, (tau[1][2] - tau0[1][2])/(tau[0][2] - tau0[0][2]), ((tau[0][2] - tau0[0][2])/l))
