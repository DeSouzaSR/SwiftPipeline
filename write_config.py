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

# Planet's names
planets_name = pd.Series(config["planets_name"])

# Semi-axis [au]
a = pd.Series([float(i) for i in config["a"]])

# Eccentricity
e = pd.Series([float(i) for i in config["e"]])

# Inclination [deg]
i = pd.Series([float(i) for i in config["i"]])

# Long. of ascending node [deg]
capom = pd.Series([float(i) for i in config["capom"]])

# arg. of peri. [deg]
omega = pd.Series([float(i) for i in config["omega"]])

# Mean anomaly [deg]
M = pd.Series([float(i) for i in config["M"]])

# Mass [kg]
mass = pd.Series([float(i) for i in config["mass"]])

# Equatorial radius [km]
radio = pd.Series([float(i) for i in config["radio"]])

# Period [days]
period = pd.Series([float(i) for i in config["period"]])


# Create a data frame from data series
planets = pd.DataFrame({'planets_name': planets_name, 'a': a, 'e': e, \
	'i': i, 'capom': capom, 'omega': omega, 'M': M, 'mass': mass, \
	'radio': radio, 'period': period})
	
# Make planets_name index
planets = planets.set_index('planets_name')

# Changing the order of the columns.
planets = planets[['a', 'e', 'i', 'capom', 'omega', 'M', 'mass', \
	'radio', 'period']]
	
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
       
simulation_name = arq.split(".")[0]
simu_suffix = glob(simulation_name + "/" + "*")
simulation_input = simulation_name + "_input"
planets.to_csv(simulation_input + "/planets.csv")

for i in simu_suffix:
    simu_clone = glob(i + "/" + "*")
    for j in simu_clone:
        with open(j + "/" + "pl.in", "w") as f:
            f.write("123\n456")
