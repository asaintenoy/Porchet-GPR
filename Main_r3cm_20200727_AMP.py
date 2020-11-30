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
import os
import pandas as pd

from picking_rada_with_amp import picked_amp
from outils import read_parameters, rada_plot
import matplotlib.pyplot as plt

#from joblib import Parallel, delayed
#os.chdir('/home/el/Codes/Porchet-GPR')
#%% Geometrie
#def des paramètres géométriques
geometry=Geometry()

#Domaine de calcul (en cm)
#Tolerance
geometry.tol=10**(-7)
# largeur
geometry.xmin=0 
geometry.xmax=40
# hauteur (elevation)
geometry.emin=0
geometry.emax = 80
# profondeur du trou en cm
geometry.dtrou = 30 # elevation du fond du trou
geometry.etrou = geometry.emax - geometry.dtrou
 # rayon du trou en cm
geometry.r=3
# hauteur d'eau imposée au fond du trou en cm
geometry.h_eau=5.0
# pas de la maille en cm
#geometry.dx = 0.1
geometry.dx = 1
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
#temps = [0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]
temps = [0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 6.00]

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

##%% Param mvg
## Teneur en eau résiduelle
#tr = [0.03]
## Teneur en eau à saturation
#ts = np.arange(0.30, 0.50, 0.05, 'float')
##ts = [0.36]
## Teneur en eau initiale
##ti = [0.07]
#ti = np.arange(0.07, 0.13, 0.03, 'float')
## Perméabilité à saturation
#lnKs = np.arange(np.log(0.005), np.log(2), 0.25, 'float')
#Ks = np.round(np.float64(np.exp(lnKs)),3)
## param fitting retention n
#n = np.arange(2, 12, 1, 'float')
## param fitting retention alpha
#alpha = np.arange(0.015, 0.085, 0.01, 'float')
#print(len(list(itertools.product(tr, ts, ti, Ks, n, alpha))))

#%% Lancement du calcul des maplitudes et écriture dans un fichier dans chaque dossier:
#%% parcours le dossier OUT
dir_name = "OUT_20200722dtrou30_rtrou3_tr5.0"

Picked = pd.DataFrame(columns=['AMP'])

for subdir in os.listdir(dir_name):
    if os.path.isfile(subdir):
        print("'%s' pas un dossier" % subdir)
        continue
    
    print(subdir)
    subdir_name = os.path.join(dir_name, subdir)
    p = read_parameters(subdir_name)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    paramMVG.porosity = paramMVG.ts
  
#    filename_TWT = os.path.join(subdir_name, 'TWT_EL.csv')
#    TWT_XP = np.genfromtxt(filename_TWT, delimiter=',',skip_header=0)            

    filename_AMP = os.path.join(subdir_name, 'AMP_EL.csv')  
    filename_radar = os.path.join(subdir_name, 'radargram__merged.out.bz2')

    if os.path.isfile(filename_radar) :
        #print("pas de radargramme zipper dans " + subdir_name)
        os.system('bunzip2 '+ filename_radar)
    
    filename_radar2 = os.path.join(subdir_name, 'radargram__merged.out')
    if not os.path.isfile(filename_radar2) :
        print('ignore ce dossier')
        continue
    
    print("Picking radargramme amplitudes dans " + subdir_name)
    
    twt_fin, amp_fin = picked_amp(filename_radar2, nT, geometry, paramMVG, paramGPRMAX, temps)   
    
    Picked['AMP'] = amp_fin
    
    # Saving Amplitude to a file
    Picked['AMP'].to_csv(filename_AMP, sep=',', encoding='utf-8', index=None, header=True)

    # zip the radargramme
    os.popen('bzip2 '+ filename_radar2)
    
    # plot des points piqués
    #figi, ax = plt.subplots()
    #ax.plot(twt_fin, amp_fin,'bo')
    #ax.plot(temps, np.abs(amp_fin[1:]),'bo')
    #ax.plot(temps, np.abs(twt_fin[1:]*10),'ro')
    #ax.grid()    
    #plt.legend(fontsize=20)
    





