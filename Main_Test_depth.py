#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 16:22:01 2020

@author: el
"""


from Forward import Forward
from param_acquisition import ParamGPRMAX, ParamMVG, Geometry
import itertools
import numpy as np
import dask
#from joblib import Parallel, delayed

#%% Geometrie
#def des paramètres géométriques
geometry=Geometry()
geometry.tol=10**(-7)
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
geometry.h_eau=5.0
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
tmax_SWMS2D = 60
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
paramGPRMAX.dx = 0.5 
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


#%% Pqrqm mvg
# Teneur en eau résiduelle
tr = 0.03
# Teneur en eau à saturation
ts = 0.36
# Teneur en eau initiale
ti = 0.09
# Perméabilité à saturation
Ks = 0.4
# param fitting retention n
n = 6
# param fitting retention alpha
alpha = 0.025
#print(len(list(itertools.product(tr, ts, ti, Ks, n, alpha))))
p=[tr,ts,ti,Ks,n,alpha]
paramMVG = ParamMVG(tr=p[0], ts=p[1],ti=p[2], Ks=p[3], n=p[4], alpha=p[5])
paramMVG.porosity = paramMVG.ts
#%% Lancement du calcul
tasks = []
dd=np.arange(25, 35, 0.5, 'float')
rr=np.arange(2,6,0.5,'float')
for p in itertools.product(dd, rr):
    # profondeur du trou en cm
    geometry.dtrou = p[0]
    geometry.r=p[1]
    # elevation du fond du trou
    geometry.etrou = geometry.emax - geometry.dtrou
    #tasks.append(dask.delayed(Forward)(geometry,paramMVG,paramGPRMAX,temps,tmax_SWMS2D))
    Forward(geometry,paramMVG,paramGPRMAX,temps,tmax_SWMS2D)

#dask.compute(tasks, scheduler='processes')



#%% test
# import time
# def connerie(ola):
#     print(str(ola))
#     time.sleep(ola)
    
    
    
# tasks = []

# for p in range(3000):#itertools.product(tr, ts, ti, Ks, n, alpha):
#     # Définition des paramètres MVG
#     #paramMVG = ParamMVG(tr=p[0], ts=p[1],ti=p[2], Ks=p[3], n=p[4], alpha=p[5])
#     #paramMVG.porosity = paramMVG.ts
#     tasks.append(dask.delayed(connerie)(p))


# dask.compute(tasks, scheduler='multiprocessing')




