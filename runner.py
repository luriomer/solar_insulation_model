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
location_name = "Las-Vegas, Nevada"
axial_tilt = 23.45*(np.pi/180) # Earth axial tilt angle [rad]
L_std = 120 # Standart longtitude line (west) [deg]
L = 115     # Actual longtitude line (west) [deg]
phi = 36*(np.pi/180) #Latitude line (north) [rad]
G0 = 1367 # Solar constant [W/m**2]
A = 0.610 # Altitude [km]
summer = [171,262] # Summer priod [start,end] [1 = Jan 1st]. In Las-Vegas June 21 - Sep 22
panel_south_angle = 32 # Panel angle to the south (about ~31 is optimal, like in Israel) [deg]
empirical_data_path = 'empirical_insulation.txt'
two_axis_tracking = False # Two axis tracking configuration. Overrides the panel angle.

''' Independent variables '''
days = np.arange(1,366,1) # Days throughout the year
hours = np.arange(0,25,1) # Hours throughout the day (already in solar time!)
delta = axial_tilt*np.sin(((days-81)*2*np.pi/365)) #Declination angle [rad]
surf_normal = surface_normal_calc(-np.tan(panel_south_angle*(np.pi/180)),0,1) # To use the degree attribute, must be yL=0 and zL=0.

''' Calling the primary engine process '''
main_calc = main_flux_calc(days,hours,phi,G0,surf_normal,delta,summer[0],summer[1],A,two_axis_tracking)
Gb = main_calc[0]
Gd = main_calc[1]
Gtot = main_calc[2]
Gav_day = main_calc[3]
Gb_av = main_calc[4]
Gd_av = main_calc[5]
Gav_hr = main_calc[6]
cos_theta = main_calc[7]
tau_b = main_calc[8]

annual = annual_calc(empirical_data_path,Gtot)
Gav_emp = annual[0]
E_sim = annual[1]
E_emp = annual[2]


''' Calling the plotter function to plot the results '''
plt.close('all')
plotter(days,hours,location_name,surf_normal,Gtot,Gav_day,Gb_av,Gd_av,Gav_emp,E_sim,E_emp,Gav_hr,two_axis_tracking,panel_south_angle)


''' Save the results to tsv .txt files inside a results folder in the local directory '''
np.savetxt("results/Gtot.txt",Gtot)