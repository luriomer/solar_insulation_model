# -*- coding: utf-8 -*-
"""
@ Omer Luria
luriomer@gmail.com 
Tel-Aviv University
Version 0.06, Dec 2017

This simulation calculates solar insulation in desired location throughout the year, and generated
graphic results.

This is the engine file.
"""
#%%
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

#%%
def surface_normal_calc(xL,yL,zL):
    ''' Surface normal direction vector. Everything in local coordinates. '''
    surf_normal_local = np.array([xL,yL,zL])
    surf_normal_local = surf_normal_local/np.linalg.norm(surf_normal_local) #Normalization
    return surf_normal_local

#%%
def Hottel_coeff(days,summer_start,summer_end,A):
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
    return [a0,a1,k]

#%%
def Klein(kT_path):
    ''' Calculation of Klein correlation according to kT datafile '''
    kT = np.loadtxt(kT_path)
    return (1.390 - 4.027*kT+5.53*kT**2-3.108*kT**3)
#%%
def solar_time_difference(days,L,Lstd):
    B = (days-1)*2*np.pi/365
    E = 229.2*(0.000075+0.001868*np.cos(B)-0.03208*np.sin(B)-0.01462*np.cos(2*B)-0.04089*np.sin(2*B))
    delta_t_solar = (4*(Lstd-L)+E)/60
    return delta_t_solar
#%%
def main_flux_calc(days,hours,phi,G0,surf_normal,delta,summer_start,summer_end,A,two_axis_tracking,Ipart,L,Lstd):
    ''' Calculation of flux in different times during the day using the correlations'''

    ''' Preperation of array variables '''
    shape = [len(days),len(hours)] # General array matrix shape [days and hours]
    sun_vec = np.zeros([3,len(days),len(hours)], dtype=float)
    omega = np.zeros(shape, dtype=float) # Hour angle
    cos_theta = np.zeros(shape, dtype=float)
    cos_zenith = np.zeros(shape, dtype=float)
    tau_b = np.zeros(shape, dtype=float)
    tau_d = np.zeros(shape, dtype=float)
    Gb = np.zeros(shape, dtype=float) # Beam flux - clear day
    Gd = np.zeros(shape, dtype=float) # Diffuse flux - clear day
    G_glob = np.zeros(shape, dtype=float) # Global flux - clear day
    qb = np.zeros(shape, dtype=float) # Beam flux - nonclear consideration
    qd = np.zeros(shape, dtype=float) # Diffuse flux - nonclear consideration
    qtot = np.zeros(shape, dtype=float) # Total flux - nonclear consideration
    qav_day = np.zeros_like(days, dtype=float) # Daily average flux
    qb_av = np.zeros_like(days, dtype=float) # Beam average flux
    qd_av = np.zeros_like(days, dtype=float) # Diffuse average flux
    qav_hr = np.zeros_like(hours, dtype=float)# Hourly average flux
    Hottel = Hottel_coeff(days,summer_start,summer_end,A)
    a0 = Hottel[0]
    a1 = Hottel[1]
    k = Hottel[2]
    
    ''' Calculation of daily parameters '''
    delta_t_solar = solar_time_difference(days,L,Lstd) #Solar noon difference [hours]
    

    ''' Main calculation loops '''
    for i in range(len(days)):
        for j in range(len(hours)):
            omega[i,j] = ((hours[j]+delta_t_solar[i]-12)*360/(24))*np.pi/180 #Hour angle [rad]
            sun_vec[0,i,j] = -np.sin(phi)*np.cos(delta[i])*np.cos(omega[i,j])+np.cos(phi)*np.sin(delta[i])
            sun_vec[1,i,j] = np.cos(delta[i])*np.sin(omega[i,j])
            sun_vec[2,i,j] = np.cos(phi)*np.cos(delta[i])*np.cos(omega[i,j])+np.sin(phi)*np.sin(delta[i])
            cos_zenith[i,j] = (np.cos(phi)*np.cos(delta[i])*np.cos(omega[i,j])+np.sin(phi)*np.sin(delta[i]))
            if two_axis_tracking:
                cos_theta[i,j] = 1.0
            else:
                cos_theta[i,j] = (np.dot(sun_vec[:,i,j],surf_normal))/(np.linalg.norm(sun_vec[:,i,j])*np.linalg.norm(surf_normal))
            if cos_zenith[i,j] > 0:
                tau_b[i,j] =max((a0[i]+a1[i]*np.exp(-k[i]/cos_zenith[i,j])),0.0)
            else:
                tau_b[i,j] = a0[i]
            Gb[i,j] = max((G0*tau_b[i,j]),0.0)
            tau_d[i,j] = max((0.271-0.294*tau_b[i,j]),0.0)
            Gd[i,j] = max((G0*tau_d[i,j]*cos_zenith[i,j]),0.0)
            G_glob[i,j] = Gb[i,j]*cos_zenith[i,j] + Gd[i,j]
            qd[i,j] = max((Ipart[i]*G_glob[i,j]),0.0)
            qb[i,j] = max((G_glob[i,j] - qd[i,j]),0.0)
            qtot[i,j] = max((qb[i,j]*cos_theta[i,j]+qd[i,j]),0.0)
        qav_day[i] = np.average(qtot[i])
        qb_av[i] = np.average(qb[i])
        qd_av[i] = np.average(qd[i])
    
    for j in range(len(hours)):
        qav_hr[j] = np.average(qtot[:,j])
        
    return (qb,qd,qtot,qav_day,qb_av,qd_av,qav_hr,cos_theta,tau_b,tau_d,cos_zenith,Ipart,delta_t_solar)
#%%

def annual_calc(empirical_path,qtot):
    ''' Importing empirical data '''
    qav_emp = np.loadtxt(empirical_path)
    ''' Calculation of total annual solar energy (average) '''
    E_sim = round(np.sum(qtot,dtype=np.float)/1e6 ,4) # Here we have all the data so simply sum it all and change it to MWh/m^2.
    E_emp = round(np.sum(qav_emp)*24/1e6,4) # Multiply every day average by 24 hours.
    return [qav_emp,E_sim,E_emp]

#%%
def plotter(days,hours,location_name,surf_normal,qtot,qav_day,qb_av,qd_av,qav_emp,E_sim,E_emp,qav_hr,two_axis_tracking,panel_south_angle):
    delta_t_noon = hours-12
    ''' Plotting the results '''
    #plt.close('all')
    
    fig1 = plt.figure()
    if two_axis_tracking:
        panel_south_angle = "Two axis tracking"
    elif panel_south_angle == 0:
        panel_south_angle = "Horizontal"
    elif panel_south_angle > 0:
        panel_south_angle = str(panel_south_angle)+"$\degree$"+" to the south"
    elif panel_south_angle < 0:
        panel_south_angle = str(abs(panel_south_angle))+"$\degree$"+" to the north"
    fig1.suptitle("Location: "+location_name+"\nSurface orientation: "
                  +str(panel_south_angle)+"\n"+"Total annual energy:"
             " E$_{simulation}$ =  "+str(E_sim)+"["+r'$\frac{MWh}{m^2}$'+"]"
             " , E$_{empirical}$ = "+str(E_emp)+"["+r'$\frac{MWh}{m^2}$'+"]" ,fontweight = "bold", fontsize=14)
    ax1 = fig1.add_subplot(221)
    ax1.plot(days,qav_day, label = "Simulation")
    ax1.plot(days,qav_emp, label = "Empirical") #Multiplting since the eimpirical averaging is for 24 hours
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Average total flux ["+r'$\frac{W}{m^2}$'+"]")
    ax1.set_title("Average daily flux (based on 24 hours)")
    ax1.legend()
    ax1.grid()
    
    ax2 = fig1.add_subplot(222)
    for i in range(int(0.5*len(delta_t_noon)-6),int(0.5*len(delta_t_noon)+6)):
        ax2.plot(days,qtot[:,i],label = "$\Delta$"+"$t_{noon}$ = "+str(delta_t_noon[i]))
        ax2.set_xlabel("Day")
        ax2.set_ylabel("Total flux ["+r'$\frac{W}{m^2}$'+"]")
        ax2.set_title("Hours through the year")
    ax2.legend(loc='right')
    ax2.grid()
    
    ax3 = fig1.add_subplot(223)
    for i in range(0,len(days),30):
        ax3.plot(hours,qtot[i,:], label ="Day ="+str(days[i]))
    ax3.set_xlabel("Hour")
    ax3.set_ylabel("Total flux ["+r'$\frac{W}{m^2}$'+"]")
    ax3.set_title("Hours through the day")
    ax3.legend()
    ax3.grid()
    
    ax4 = fig1.add_subplot(224)
    ax4.plot(days,qb_av, label = "Beam")
    ax4.plot(days,qd_av, label = "Diffuse")
    ax4.set_xlabel("Day")
    ax4.set_ylabel("Average flux components ["+r'$\frac{W}{m^2}$'+"]")
    ax4.set_title("Average beam vs. diffuse flux")
    ax4.legend()
    ax4.grid()
    return