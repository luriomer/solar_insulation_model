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
axial_tilt = 23.45*np.pi/180 # Earth axial tilt angle [rad]
L_std = 120 # Standart longtitude line (west) [deg]
L = 115     # Actual longtitude line (west) [deg]
phi = 36*np.pi/180 #Latitude line (north) [rad]
G0 = 1367 # Solar constant [W/m**2]
A = 0.610 # Altitude [km]

''' Independent variables '''
N = np.linspace(1,365,365) # Days throughout the year
B = (2*np.pi*(N-1)/365)*np.pi/180 # Assisting parameter [rad]
E = 229.2*(0.000075+0.001868*np.cos(B)-0.03208*np.sin(B)-0.01462*np.cos(2*B)-0.04089*np.sin(2*B))

''' Calculation of declination angle '''
delta = axial_tilt*np.sin((B-80*2*np.pi/365)*np.pi/180) #Declination angle [rad]
#delta_t_solar = 4*(L_std-L)+E #Solar time difference (noon to noon) [minutes]

''' Calculation of Hottel correlation coefficients. 
    In Las Vegas, summer is 21 June to 22 September '''
r0 = np.zeros_like(N)
r1 = np.zeros_like(N)
rk = np.zeros_like(N)
for i in range(len(N)):
    if (N[i]>=171 and N[i]<=262):
        r0[i] = 1.03
        r1[i] = 1.01
        rk[i] = 1.00
    else:
        r0[i] = 0.97
        r1[i] = 0.99
        rk[i] = 1.02
a0 = r0*(0.4237-0.00821*(6-A)**2)
a1 = r1*(0.5055+0.00595*(6.5-A)**2)
k = rk*(0.2711+0.01858*(2.5-A)**2)

''' Calculation and plot of flux in different times during the day using the correlations'''
Gan = zeros_like(N)
plt.close('all')
fig, ax1 = plt.subplots(2, sharex = True)
fig, ax2 = plt.subplots()
delta_ts = np.arange(-5,6,1)
for delta_t_solar in delta_ts:
    omega = delta_t_solar*360/(24*60) #Hour angle [deg]
    cos_theta_z = (np.cos(phi)*np.cos(delta)*np.cos(omega)+np.sin(phi)*np.sin(delta))
    AM = 1/cos_theta_z
    tau_b = a0+a1*np.exp(-k*AM)
    Gb = G0*tau_b
    tau_d = cos_theta_z*(0.271-0.294*tau_b)
    Gd = G0*tau_d
    Gtot = Gb + Gd
    ax1[0].plot(N,Gb, label = 'Solar noon difference = '+str(np.abs(delta_t_solar)))
    ax1[1].plot(N,Gd, label = 'Solar noon difference = '+str(np.abs(delta_t_solar)))
    ax2.plot(N,Gtot, label = 'Solar noon difference = '+str(np.abs(delta_t_solar)))
#ax1[0].set_xlabel('Day number N (N=1 for 1 January)')
ax1[0].set_ylabel('Beam insulation flux G [W/m**2]')
ax1[0].legend(loc = 4)
ax1[0].grid()
ax1[1].set_xlabel('Day number N (N=1 for 1 January)')
ax1[1].set_ylabel('Diffuse insulation flux G [W/m**2]')
ax1[1].legend(loc = 4)
ax1[1].grid()
ax2.set_xlabel('Day number N (N=1 for 1 January)')
ax2.set_ylabel('Total insulation flux G [W/m**2]')
ax2.legend(loc = 4)
ax2.grid()