# -*- coding: utf-8 -*-
"""
@ Omer Luria
luriomer@gmail.com 
Tel-Aviv University
Version 0.05, Dec 2017

This simulation calculates solar insulation in desired location throughout the year, and generated
graphic results.

This is the runner file. Run it with the desired parameters, in the same folder with the engine.
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from model_engine import surface_normal_calc,main_flux_calc,annual_calc,plotter


''' Analysis parameters: change as neccessary'''
location_name = "Las-Vegas, NV"
axial_tilt = 23.45*(np.pi/180) # Earth axial tilt angle [rad]
L_std = 120 # Standart longtitude line (west) [deg]
L = 115     # Actual longtitude line (west) [deg]
phi = 36*(np.pi/180) #Latitude line (north) [rad]
G0 = 1367 # Solar constant [W/m**2]
A = 0.610 # Altitude [km]
summer_start = 171 # First day of summer [days, 1 = Jan 1st]. In Las-Vegas June 21
summer_end = 262 # Last day of summer [days, 1 = Jan 1st]. In Las-Vegas Sep 22

''' Independent variables '''
days = np.arange(1,366,1) # Days throughout the year
hours = np.arange(7,18,1) # Daylight hours
B = (2*np.pi*(days-1)/365)*np.pi/180 # Assisting parameter [rad]

''' Calling surface normal calculation function '''
surf_normal = surface_normal_calc(0,2,0)


''' Calculation of declination angle '''
delta = axial_tilt*np.sin((B-80*2*np.pi/365)*np.pi/180) #Declination angle [rad]

''' Calling the primary engine process '''
main_calc = main_flux_calc(days,hours,phi,G0,surf_normal,delta,summer_start,summer_end,A)#,a0,a1,k)
Gb = main_calc[0]
Gd = main_calc[1]
Gtot = main_calc[2]
Gav_day = main_calc[3]
Gb_av = main_calc[4]
Gd_av = main_calc[5]
Gav_hr = main_calc[6]

annual = annual_calc('empirical_insulation.txt',Gtot)
Gav_emp = annual[0]
E_sim = annual[1]
E_emp = annual[2]

''' Calling the plotter function to plot the results '''
plt.close('all')
plotter(days,hours,location_name,surf_normal,Gtot,Gav_day,Gb_av,Gd_av,Gav_emp,E_sim,E_emp,Gav_hr)