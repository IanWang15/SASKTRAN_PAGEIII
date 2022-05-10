import sasktran as sk
from sasktran.geometry import VerticalImage
import numpy as np
import matplotlib.pyplot as plt

# altitude follows OMPS
altarr = np.arange(0.5,40.6,1.)
#print(altarr)

#import pdb
#pdb.set_trace()
#print('stop here')

# First recreate our geometry and atmosphere classes
geometry = VerticalImage()
geometry.from_sza_saa(sza=30, saa=60, lat=0, lon=0, tanalts_km=altarr, mjd=54372, locallook=0,
                      satalt_km=840)#, refalt_km=15)

sizelist = [0.1,0.2,0.5,1.]
nclist = [0.01,0.1,1,10,100]

rows = []

for size in sizelist:
    for inc in nclist:
        print(size,inc)
        atmosphere = sk.Atmosphere()
# Turn off Rayleigh scattering and gaseous in the air
#atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
#atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        atmosphere.atmospheric_state = sk.MSIS90()
#atmosphere['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
#atmosphere['no2'] = sk.Species(sk.NO2Vandaele1998(), sk.Pratmo())
#configure_for_cloud(grid_spacing_m: float, cloud_altitudes: numpy.array, cloud_diffuse_spacing_m: float = 100, max_optical_depth_of_cell: float = 0.01, min_extinction_ratio_of_cell: float = 1)

        particle_size_dist = sk.ClimatologyUserDefined(altitudes=[0, 16000], \
              values={'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': [1.6, 1.6], \
              'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': [size, size]})

        aerosol_optprop = sk.MieAerosol(particlesize_climatology=particle_size_dist, species='H2SO4')

        xs = aerosol_optprop.calculate_cross_sections(sk.MSIS90(),latitude=0, longitude=0, altitude=10000, mjd=54372, wavelengths=[510.,600.,675.,745.,868.,997.]).total
        print('xs ',xs)

        xs = xs * inc

        print(np.shape(xs))
        print(xs[:])
        row = xs[:]
        row = np.insert(row, 0, inc)
        row = np.insert(row, 0, size)
        rows.append(row)
        print(row)
        print('-------------')
        print(rows)

import csv

# field names
fields = ['Size', 'N', 510., 600., 675., 745., 868., 997.]

# name of csv file
filename = "tmpxsaer.csv"
# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)

    # writing the fields
    csvwriter.writerow(fields)

    # writing the data rows
    csvwriter.writerows(rows)

