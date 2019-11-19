#!/usr/bin/env python3
from netCDF4 import Dataset # pip install netCDF4

in_file = 'spifs.nc'
rootgrp = Dataset(in_file, "r")

# Separate superparameterized grid points from non-superparameterized ones.
# Superparameterized grid points have more variables.
SP_groups = []
extra_groups = []
for g in rootgrp.groups.keys():
    if len(rootgrp[g].variables) > 40:
        SP_groups.append(g)
    else:
        extra_groups.append(g)
print ('SP groups:', SP_groups)

# Print profiles from one group
group = SP_groups[0]
time_index = 5
lat = rootgrp[group+'/lat'][:]
lon = rootgrp[group+'/lon'][:]
time = rootgrp['/Time'][time_index]
time_unit = rootgrp['/Time'].units
print ('group:', group, 'lat:', lat, 'lon:', lon, 'at time', time, time_unit)
print ()
Zf = rootgrp[group+'/Zf'][time_index][:]   # Height
T  = rootgrp[group+'/T'][time_index][:]    # Temperature
SH  = rootgrp[group+'/SH'][time_index][:]  # Specific humidity
U  = rootgrp[group+'/U'][time_index][:]    # East-West wind
V  = rootgrp[group+'/V'][time_index][:]    # North-South wind

for i in range(len(Zf)):
    print ('%8.1f %5.1f %6.4f %4.1f %4.1f'%(Zf[i], T[i], SH[i], U[i], V[i]))



