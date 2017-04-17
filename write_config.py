#!/usr/bin/env python

from glob import glob
from configobj import ConfigObj
import numpy as np
import pandas as pd
import math
import oe2pv


# Create planets data (J2000)
# 
# Data from https://nssdc.gsfc.nasa.gov/planetary/
# Data sequence: 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus',
#  'Neptune'.

arq = "Mercury.ini" #sys.argv[1]
config = ConfigObj(arq)

simulation_name = arq.split(".")[0]
simulation_input = simulation_name + "_input"
planets = pd.read_csv(simulation_input + "/planets_ini.csv", index_col = "planets_name")

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
    x[j], y[j], z[j], vx[j], vy[j], vz[j] = oe2pv.orbel_el2xv(gm[j],\
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
       

simu_suffix = glob(simulation_name + "/" + "*")
planets.to_csv(simulation_input + "/planets_input.csv")

for i in simu_suffix:
    simu_clone = glob(i + "/" + "*")
    for j in simu_clone:
        with open(j + "/" + "pl.in", "w") as f:
            f.write("123\n456\n789")
