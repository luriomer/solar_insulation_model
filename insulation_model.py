# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 17:39:30 2017

@author: Omer Luria
Version 0.03

This simulation calculates solar insulation in desired location throughout the year, and generated
graphic results. In this case the location chosen was Las Vegas, NV.
The current model takes into account constant surface direction.
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


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


''' Surface normal direction vector. Everything in local coordinates. '''
surf_normal = np.array([0,1,0]) # Change according to the definition of coordinates.
surf_normal = surf_normal/np.sqrt(np.dot(surf_normal,surf_normal)) #Normalization of surface vector


''' Independent variables '''
days = np.arange(1,366,1) # Days throughout the year
hours = np.arange(7,18,1)
B = (2*np.pi*(days-1)/365)*np.pi/180 # Assisting parameter [rad]


''' Calculation of declination angle '''
delta = axial_tilt*np.sin((B-80*2*np.pi/365)*np.pi/180) #Declination angle [rad]


''' Calculation of Hottel correlation coefficients '''
r0 = np.zeros_like(days, dtype=float)
r1 = np.zeros_like(days, dtype=float)
rk = np.zeros_like(days, dtype=float)
for i in range(len(days)):
    if (days[i]>=summer_start and days[i]<=summer_end):
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
sun_vec = np.zeros([3,len(days),len(hours)])
cos_theta = np.zeros(shape, dtype=float)
tau_b = np.zeros(shape, dtype=float)
tau_d = np.zeros(shape, dtype=float)
Gb = np.zeros(shape, dtype=float) # Beam flux (instant)
Gd = np.zeros(shape, dtype=float) # Diffuse flux (instant)
Gtot = np.zeros(shape, dtype=float) # Total flux (instant)
Gav_day = np.zeros_like(days, dtype=float) # Daily average flux
Gb_av = np.zeros_like(days, dtype=float) # Beam average flux
Gd_av = np.zeros_like(days, dtype=float) # Diffuse average flux
Gav_hr = np.zeros_like(hours, dtype=float)# Hourly average flux


''' Calculation of daily parameters '''
delta_t_solar = hours-12
omega = (delta_t_solar*360/(24))*np.pi/180 #Hour angle [rad]


''' Calculation of flux in different times during the day using the correlations'''

for i in range(len(days)):
    for j in range(len(delta_t_solar)):
        sun_vec[0,i,j] = -np.sin(phi)*np.cos(delta[i])*np.cos(omega[j])+np.cos(phi)*np.sin(delta[i])
        sun_vec[1,i,j] = np.cos(delta[i])*np.sin(omega[j])
        sun_vec[2,i,j] = np.cos(phi)*np.cos(delta[i])*np.cos(omega[j])+np.sin(phi)*np.sin(delta[i])
        cos_theta[i,j] = (np.dot(sun_vec[:,i,j],surf_normal))/(np.linalg.norm(sun_vec[:,i,j])*np.linalg.norm(surf_normal))
        #zenith[i,j] = (np.cos(phi)*np.cos(delta[i])*np.cos(omega[j])+np.sin(phi)*np.sin(delta[i]))
        tau_b[i,j] =max((a0[i]+a1[i]*np.exp(-k[i]/cos_theta[i,j])),0)
        Gb[i,j] = max((G0*tau_b[i,j]*cos_theta[i,j]),0)
        tau_d[i,j] = max((cos_theta[i,j]*(0.271-0.294*tau_b[i,j])),0)
        Gd[i,j] = max((G0*tau_d[i,j]*cos_theta[i,j]),0)
        Gtot[i,j] = Gb[i,j] + Gd[i,j]
    Gav_day[i] = np.average(Gtot[i])
    Gb_av[i] = np.average(Gb[i])
    Gd_av[i] = np.average(Gd[i])


''' Importing empirical data '''
Gav_emp = np.loadtxt('empirical_insulation.txt')


''' Calculation of total annual solar energy (average) '''
E_sim = round(np.sum(Gtot)/1e6 ,2) # Here we have all the data so simply sum it all and change it to MWh.
E_emp = round(np.sum(Gav_emp)*24/1e6,2) # Multiply every day average by 12 hours of insulation.


''' Plotting the results '''
plt.close('all')

fig1 = plt.figure()
fig1.suptitle("Simulation results for "+location_name+"\nSurface normal (local coordinates) = "
              +str(surf_normal), fontweight = "bold", fontsize=18)
ax1 = fig1.add_subplot(221)
ax1.plot(days,Gav_day, label = "Simulation")
ax1.plot(days,Gav_emp*2, label = "Empirical") #Multiplting since the eimpirical averaging is for 24 hours
ax1.set_xlabel("Day")
ax1.set_ylabel("Average flux [W/m$^2$]")
ax1.set_title("Average daily flux based on 12 hours of daylight")
ax1.annotate("Total annual energy per unit area: \n "+
         "--------------------------------------------------\n"+
         "Theoretical (simulation) E = "+str(E_sim)+" [MWh/m$^2$] \n"+
         "Empirical (atmospheric data) E = "+str(E_emp)+" [MWh/m$^2$]", 
         fontsize = 11, color="red", fontweight = "bold", xy=(70, 240))
ax1.legend()
ax1.grid()

ax2 = fig1.add_subplot(222)
for i in range(len(delta_t_solar)):
    ax2.plot(days,Gtot[:,i],label = "$\Delta$t solar = "+str(delta_t_solar[i]))
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Total flux [W/m$^2$]")
    ax2.set_title("Total daily flux")
    ax2.legend()
    ax2.grid()

ax3 = fig1.add_subplot(223)
for j in range(len(hours)):
    Gav_hr[j] = np.average(Gtot[:,j])
ax3.plot(hours,Gav_hr)
ax3.set_xlabel("Hour")
ax3.set_ylabel("Total flux [W/m$^2$]")
ax3.set_title("Average hourly flux")
ax3.legend()
ax3.grid()

ax4 = fig1.add_subplot(224)
ax4.plot(days,Gb_av, label = "Beam")
ax4.plot(days,Gd_av, label = "Diffuse")
ax4.set_xlabel("Day")
ax4.set_ylabel("Flux components [W/m$^2$]")
ax4.set_title("Beam vs. Diffuse flux")
ax4.legend()
ax4.grid()
