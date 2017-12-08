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
#%%
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

#%%
def surface_normal_calc(xL,yL,zL):
    ''' Surface normal direction vector. Everything in local coordinates. '''
    surf_normal_local = np.array([xL,yL,zL])/np.sqrt(xL**2+yL**2+zL**2) #Normalization
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
def main_flux_calc(days,hours,phi,G0,surf_normal,delta,summer_start,summer_end,A,two_axis_tracking):
    ''' Calculation of flux in different times during the day using the correlations'''

    ''' Preperation of array variables '''
    shape = [len(days),len(hours)] # General array matrix shape [days and hours]
    sun_vec = np.zeros([3,len(days),len(hours)])
    cos_theta = np.zeros(shape, dtype=float)
    cos_zenith = np.zeros(shape, dtype=float)
    tau_b = np.zeros(shape, dtype=float)
    tau_d = np.zeros(shape, dtype=float)
    Gb = np.zeros(shape, dtype=float) # Beam flux (instant)
    Gd = np.zeros(shape, dtype=float) # Diffuse flux (instant)
    Gtot = np.zeros(shape, dtype=float) # Total flux (instant)
    Gav_day = np.zeros_like(days, dtype=float) # Daily average flux
    Gb_av = np.zeros_like(days, dtype=float) # Beam average flux
    Gd_av = np.zeros_like(days, dtype=float) # Diffuse average flux
    Gav_hr = np.zeros_like(hours, dtype=float)# Hourly average flux
    Hottel = Hottel_coeff(days,summer_start,summer_end,A)
    a0 = Hottel[0]
    a1 = Hottel[1]
    k = Hottel[2]
    
    ''' Calculation of daily parameters '''
    delta_t_solar = hours-12
    omega = (delta_t_solar*360/(24))*np.pi/180 #Hour angle [rad]
    

    ''' Main calculation loops '''
    for i in range(len(days)):
        for j in range(len(delta_t_solar)):
            sun_vec[0,i,j] = -np.sin(phi)*np.cos(delta[i])*np.cos(omega[j])+np.cos(phi)*np.sin(delta[i])
            sun_vec[1,i,j] = np.cos(delta[i])*np.sin(omega[j])
            sun_vec[2,i,j] = np.cos(phi)*np.cos(delta[i])*np.cos(omega[j])+np.sin(phi)*np.sin(delta[i])
            cos_zenith[i,j] = (np.cos(phi)*np.cos(delta[i])*np.cos(omega[j])+np.sin(phi)*np.sin(delta[i]))
            if two_axis_tracking:
                cos_theta[i,j] = cos_zenith[i,j]
            else:
                cos_theta[i,j] = (np.dot(sun_vec[:,i,j],surf_normal))/(np.linalg.norm(sun_vec[:,i,j])*np.linalg.norm(surf_normal))
            tau_b[i,j] =max((a0[i]+a1[i]*np.exp(-k[i]/cos_zenith[i,j])),0)
            Gb[i,j] = max((G0*tau_b[i,j]*cos_theta[i,j]),0)
            tau_d[i,j] = max((cos_zenith[i,j]*(0.271-0.294*tau_b[i,j])),0)
            Gd[i,j] = max((G0*tau_d[i,j]*cos_theta[i,j]),0)
            Gtot[i,j] = Gb[i,j] + Gd[i,j]
        Gav_day[i] = np.average(Gtot[i,:])
        Gb_av[i] = np.average(Gb[i])
        Gd_av[i] = np.average(Gd[i])
    
    for j in range(len(hours)):
        Gav_hr[j] = np.average(Gtot[:,j])
        
    return (Gb,Gd,Gtot,Gav_day,Gb_av,Gd_av,Gav_hr)
#%%

def annual_calc(empirical_path,Gtot):
    ''' Importing empirical data '''
    Gav_emp = np.loadtxt(empirical_path)
    ''' Calculation of total annual solar energy (average) '''
    E_sim = round(np.sum(Gtot,dtype=np.float)/1e6 ,2) # Here we have all the data so simply sum it all and change it to MWh/m^2.
    E_emp = round(np.sum(Gav_emp)*24/1e6,2) # Multiply every day average by 24 hours.
    return [Gav_emp,E_sim,E_emp]

#%%
def plotter(days,hours,location_name,surf_normal,Gtot,Gav_day,Gb_av,Gd_av,Gav_emp,E_sim,E_emp,Gav_hr,two_axis_tracking):
    delta_t_solar=hours-12
    ''' Plotting the results '''
    #plt.close('all')
    
    fig1 = plt.figure()
    if two_axis_tracking:
        surf_normal = "Two axis surface"
    fig1.suptitle("Location: "+location_name+"\nSurface normal (local coordinates) --> "
                  +str(surf_normal)+"\n"+"Total annual energy:"
             " E$_{simulation}$ =  "+str(E_sim)+"["+r'$\frac{MWh}{m^2}$'+"]"
             " , E$_{empirical}$ = "+str(E_emp)+"["+r'$\frac{MWh}{m^2}$'+"]" ,fontweight = "bold", fontsize=12)
    ax1 = fig1.add_subplot(221)
    ax1.plot(days,Gav_day, label = "Simulation")
    ax1.plot(days,Gav_emp, label = "Empirical") #Multiplting since the eimpirical averaging is for 24 hours
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Average global flux ["+r'$\frac{W}{m^2}$'+"]")
    ax1.set_title("Average daily flux (based on 24 hours)")
    ax1.legend()
    ax1.grid()
    
    ax2 = fig1.add_subplot(222)
    for i in range(int(0.5*len(delta_t_solar)-6),int(0.5*len(delta_t_solar)+6)):
        ax2.plot(days,Gtot[:,i],label = "$\Delta$"+"$t_{solar}$ = "+str(delta_t_solar[i]))
        ax2.set_xlabel("Day")
        ax2.set_ylabel("Global flux ["+r'$\frac{W}{m^2}$'+"]")
        ax2.set_title("Hours through the year")
        ax2.legend()
        ax2.grid()
    
    ax3 = fig1.add_subplot(223)
    for i in range(0,len(days),30):
        ax3.plot(hours,Gtot[i,:], label ="Day ="+str(days[i]))
    ax3.set_xlabel("Solar Hour")
    ax3.set_ylabel("Global flux ["+r'$\frac{W}{m^2}$'+"]")
    ax3.set_title("Hours through the day")
    ax3.legend()
    ax3.grid()
    
    ax4 = fig1.add_subplot(224)
    ax4.plot(days,Gb_av, label = "Beam")
    ax4.plot(days,Gd_av, label = "Diffuse")
    ax4.set_xlabel("Day")
    ax4.set_ylabel("Average flux components ["+r'$\frac{W}{m^2}$'+"]")
    ax4.set_title("Average beam vs. diffuse flux")
    ax4.legend()
    ax4.grid()
    return
#%%