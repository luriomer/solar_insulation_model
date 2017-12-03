# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 17:39:30 2017

@author: Omer Luria
Version 0.03

This simulation calculates solar insulation in desired location throughout the year, and generated
graphic results. In this case the location chosen was Las Vegas, NV.
The current model takes into account constant surface direction.
"""
#%%
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from model_engine import surface_normal_calc,main_flux_calc,annual_calc,plotter

#%%
''' Las Vegas and universal parameters '''
location_name = "Las-Vegas, NV"
axial_tilt = 23.45*(np.pi/180) # Earth axial tilt angle [rad]
L_std = 120 # Standart longtitude line (west) [deg]
L = 115     # Actual longtitude line (west) [deg]
phi = 36*(np.pi/180) #Latitude line (north) [rad]
G0 = 1367 # Solar constant [W/m**2]
A = 0.610 # Altitude [km]
summer_start = 171 # First day of summer [days, 1 = Jan 1st]. In Las-Vegas June 21
summer_end = 262 # Last day of summer [days, 1 = Jan 1st]. In Las-Vegas Sep 22

#%%
surf_normal = surface_normal_calc(0,2,0)

#%%
''' Independent variables '''
days = np.arange(1,366,1) # Days throughout the year
hours = np.arange(7,18,1) # Daylight hours
B = (2*np.pi*(days-1)/365)*np.pi/180 # Assisting parameter [rad]

#%%
''' Calculation of declination angle '''
delta = axial_tilt*np.sin((B-80*2*np.pi/365)*np.pi/180) #Declination angle [rad]

#%%
main_calc = main_flux_calc(days,hours,phi,G0,surf_normal,delta,summer_start,summer_end,A)#,a0,a1,k)
Gb = main_calc[0]
Gd = main_calc[1]
Gtot = main_calc[2]
Gav_day = main_calc[3]
Gb_av = main_calc[4]
Gd_av = main_calc[5]
Gav_hr = main_calc[6]

#%%
annual = annual_calc('empirical_insulation.txt',Gtot)
Gav_emp = annual[0]
E_sim = annual[1]
E_emp = annual[2]

#%%
plt.close('all')
plotter(days,hours,location_name,surf_normal,Gtot,Gav_day,Gb_av,Gd_av,Gav_emp,E_sim,E_emp,Gav_hr)