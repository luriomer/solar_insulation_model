# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 17:39:30 2017

@author: Omer Luria
Version 0.01

This simulation calculates solar insulation in desired location throughout the year, and generated
graphic results. In this case the location chosen was Las Vegas, NV.
The results currently take into account horizontal surface alone.

"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


''' Las Vegas and universal parameters '''
axial_tilt = 23.45*(np.pi/180) # Earth axial tilt angle [rad]
L_std = 120 # Standart longtitude line (west) [deg]
L = 115     # Actual longtitude line (west) [deg]
phi = 36*(np.pi/180) #Latitude line (north) [rad]
G0 = 1367 # Solar constant [W/m**2]
A = 0.610 # Altitude [km]


''' Independent variables '''
days = np.arange(1,366,1) # Days throughout the year
hours = np.arange(7,18,1)
B = (2*np.pi*(days-1)/365)*np.pi/180 # Assisting parameter [rad]


''' Calculation of declination angle '''
delta = axial_tilt*np.sin((B-80*2*np.pi/365)*np.pi/180) #Declination angle [rad]


''' Calculation of Hottel correlation coefficients. 
    In Las Vegas, summer is 21 June to 22 September '''
r0 = np.zeros_like(days, dtype=float)
r1 = np.zeros_like(days, dtype=float)
rk = np.zeros_like(days, dtype=float)
for i in range(len(days)):
    if (days[i]>=171 and days[i]<=262):
        r0[i] = 0.97
        r1[i] = 0.99
        rk[i] = 1.02
    else:
        r0[i] = 1.03
        r1[i] = 1.01
        rk[i] = 1.00
a0 = r0*(0.4237-0.00821*(6-A)**2)
a1 = r1*(0.5055+0.00595*(6.5-A)**2)
k = rk*(0.2711+0.01858*(2.5-A)**2)


''' Preperation of array variables '''
shape = [len(days),len(hours)] # General array matrix shape [days and hours]
cos_theta_z = np.zeros(shape, dtype=float)
tau_b = np.zeros(shape, dtype=float)
tau_d = np.zeros(shape, dtype=float)
Gb = np.zeros(shape, dtype=float) # Beam flux (instant)
Gd = np.zeros(shape, dtype=float) # Diffuse flux (instant)
Gtot = np.zeros(shape, dtype=float) # Total flux (instant)
Gav = np.zeros_like(days, dtype=float) # Daily total flux (average)


''' Calculation of daily parameters '''
delta_t_solar = hours-12
omega = (delta_t_solar*360/(24))*np.pi/180 #Hour angle [rad]


''' Calculation of flux in different times during the day using the correlations'''

for i in range(len(days)):
    for j in range(len(delta_t_solar)):
        cos_theta_z[i,j] = (np.cos(phi)*np.cos(delta[i])*np.cos(omega[j])+np.sin(phi)*np.sin(delta[i]))
        tau_b[i,j] = a0[i]+a1[i]*np.exp(-k[i]/cos_theta_z[i,j])
        Gb[i,j] = G0*tau_b[i,j]
        tau_d[i,j] = cos_theta_z[i,j]*(0.271-0.294*tau_b[i,j])
        Gd[i,j] = G0*tau_d[i,j]
        Gtot[i,j] = Gb[i,j] + Gd[i,j]
    Gav[i] = np.average(Gtot[i])


''' Importing empirical data '''
Gav_emp = np.loadtxt('empirical_insulation.txt')


''' Plotting the results '''
plt.close('all')

fig1 = plt.figure()

ax1 = fig1.add_subplot(221)
ax1.plot(days,Gav, label = "Simulation")
ax1.plot(days,Gav_emp, label = "Empirical")
ax1.set_xlabel("Day")
ax1.set_ylabel("Average flux [W/m$^2$]")
ax1.set_title("Average daily throught the year")
ax1.legend()
ax1.grid()

ax2 = fig1.add_subplot(222)
for i in range(len(delta_t_solar)):
    ax2.plot(days,Gtot[:,i],label = "$\Delta$t = "+str(delta_t_solar[i]))
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Total flux [W/m$^2$]")
    ax2.set_title("Total daily flux throughout the year")
    ax2.legend()
    ax2.grid()

ax3 = fig1.add_subplot(223)
for j in range(len(days)):
    ax3.plot(hours,Gtot[j])
    ax3.set_xlabel("Hour")
    ax3.set_ylabel("Total flux [W/m$^2$]")
    ax3.set_title("Total flux throughout the day")
    ax3.legend()
    ax3.grid()