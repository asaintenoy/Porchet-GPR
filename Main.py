#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:22:11 2020

@author: el
"""
import matplotlib.pyplot as plt
from param_acquisition import Geometry,ParamMVG, ParamGPRMAX
import os
os.chdir('/home/el/Codes/Porchet-GPR')

#%% Param MVG
# Teneur en eau résiduelle
tr = 0.0105
# Teneur en eau à saturation
ts = 0.35
# Teneur en eau initiale
ti = 0.07
# Perméabilité à saturation
Ks = 0.215
# param fitting retention n
n = 5
# param fitting retention alpha
alpha = 0.03
pVg=ParamMVG(tr=tr, ts=ts, ti=ti, Ks=Ks, n=n, alpha=alpha)
pVg.porosity = pVg.ts

#%% Geometrie
#def des paramètres géométriques
geometry=Geometry()

#Domaine de calcul (en cm)
# largeur
geometry.xmin=0 
geometry.xmax=40
# hauteur (elevation)
geometry.emin=0
geometry.emax = 80
# profondeur du trou en cm
geometry.dtrou = 30
# elevation du fond du trou
geometry.etrou = geometry.emax - geometry.dtrou
 # rayon du trou en cm
geometry.r=4
# hauteur d'eau imposée au fond du trou en cm
geometry.h_eau=10#5.0
# pas de la maille en cm
geometry.dx = 0.1
#geometry.dx = 1
# profondeur sous le trou (cm) jusqu'où on souhaite un maillage affiné. 
geometry.zaff= 20
#largeur horizontal de la zone affinée (cm)
geometry.waff=20
# elevation de l'affinage
geometry.eaff=geometry.etrou-geometry.zaff 
# contrainte d'angle min pour mesh 
geometry.quality=33
# maximum triangle size  (m*²)
geometry.area=5
# tupple for mesh generation 
geometry.smooth=[1,5]

#%% def de paramètres SWMS
#Temps d'infiltration où à lieu le calcul de chaque trace  (minutes) 
#temps=[1.00, 2.00]
temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

# Temps max de calcul SWMS2D au delà duquel on arrète le calcul (secondes)
tmax_SWMS2D = 600
#tmax_SWMS2D = 10

nT=len(temps)

#%% Definition des param gprMax
paramGPRMAX=ParamGPRMAX()
# Domaine de calcul (cm)
paramGPRMAX.xmin = geometry.xmin
paramGPRMAX.xmax = geometry.xmax
paramGPRMAX.zmin = geometry.emin
paramGPRMAX.zmax = geometry.emax
# Taille des mailles (cm)
paramGPRMAX.dx = 1.0 
# Electrical conductivity of the medium
paramGPRMAX.sigma=0.0000
# Relative dielectric permittivity of water
paramGPRMAX.eps_w=80.1
# Relative dielectric permittivity of PVC
paramGPRMAX.eps_pvc=3
# Relative dielectric permittivity of pure silice
paramGPRMAX.eps_s=2.5
# Ricker signal central frequency (Hz)
paramGPRMAX.wave_freq = 1000e6
# Frequence max du signal EM (Hz)
paramGPRMAX.freq_max = 2.8 * paramGPRMAX.wave_freq
# Distance between hole middle and source (m)
paramGPRMAX.d_emet = 0.18
# Distance between hole middle and receiving antenna (m)
paramGPRMAX.d_recept = 0.22
# param qui raffine le pas spatial (par défaut 10 d'après doc gprmax)
paramGPRMAX.spatial_step = 5
# Trace time window (ns)
paramGPRMAX.time = 30e-9
#time_step_stability_factor (pas utilisé pour le moment...)
paramGPRMAX.fac_dt = 0.2 
#%% Forward
from Forward import Forward
[TWT,Vol]=Forward(geometry,pVg,paramGPRMAX,temps,tmax_SWMS2D)

#import multiprocessing as mp
#print("Number of processors: ", mp.cpu_count())

#%%
#temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

XX=temps.copy()
XX.insert(0, 0)
plt.close('all')
fig,ax1=plt.subplots(1,1,figsize=(25,15))
ax1.plot(XX,TWT,'r',label='TWT(ns)')
ax1.grid()
ax1.set_xlabel('Minutes')
ax1.set_ylabel('TWT(ns)')
ax2 = ax1.twinx()
ax2.plot(XX,Vol,'b',label='Volume')
ax2.grid()
ax2.set_ylabel('Volume(mL)')
fig.legend()

#%%
from outils import rada_plot
import h5py
filepath='OUTTESRdtrou30_rtrou2.5_tr10/ts0.35_ti0.07_tr0.0105_n5_alpha0.03_Ks0.215/'
filename='radargram__merged.out'
f = h5py.File(filepath + filename, 'r')

path = '/rxs/rx1/'
data = f['%s%s' % (path, 'Ez')][:,:]
plt.figure(figsize=(10,15))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.imshow(data[:,:],aspect=0.005)