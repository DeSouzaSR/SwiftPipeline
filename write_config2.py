#!/usr/bin/env python
#-*- coding:utf-8 -*-

import numpy as np
import pandas as pd
import math
from oe2pv import orbel_el2xv


planets = pd.read_csv("Mercury_input/planets_ini.csv", index_col = "planets_name")

#Create new column, considering G = 1
# Mass of the Sum, in kg
mass_sun_kg = 1988500e24

# Mass of the Sun, with G = 1
mass_sun_grav = 2.959139768995959e-04

# Conic section is ellipse
ialpha = -1

# Gravitational factor of the Sun
gm =  2.959139768995959e-04

# Create mass_grav column
planets['mass_grav'] = planets.mass * mass_sun_grav / mass_sun_kg

# Create gmpl
planets['gmpl'] = gm + planets.mass_grav

# Creating variables to use the "orbel" function
gm = planets['gmpl']
a = planets['a']
e = planets['e']
inc = planets['i']
capom = planets['capom']
omega = planets['omega']
capm = planets['M']
P = planets['period']
rpl = planets['radio']

# Create position and velocity columns
#
# The module eo2pv.so, constructed from the Swift conversion subroutine,
# will be used.

len_planets = len(planets)

x = np.zeros(len_planets)
y = np.zeros(len_planets)
z = np.zeros(len_planets)
vx = np.zeros(len_planets)
vy = np.zeros(len_planets)
vz = np.zeros(len_planets)

for j in range(len(planets)):
    x[j], y[j], z[j], vx[j], vy[j], vz[j] = orbel_el2xv(gm[j],\
										ialpha, a[j],e[j],\
        									math.radians(inc[j]),\
                                            math.radians(capom[j]),\
                                            math.radians(omega[j]),\
                                            capm[j])


# Create colums x, y, v, vx, vy and vz
planets['x'] = x
planets['y'] = y
planets['z'] = z
planets['vx'] = vx
planets['vy'] = vy
planets['vz'] = vz

print(planets)