#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:22:52 2020

@author: sainteno
"""
#%%
import math
from collections import namedtuple
from maillage_GPRMAX import CRIM

#Def des paramètres de géométrie du modéle
class Geometry :
    def __repr__(self):
        return "dtrou{dtrou}_rtrou{r}_tr{h_eau}".format(dtrou=self.dtrou,r=self.r,h_eau=self.h_eau)        

# Definition des paramètres MVG
class ParamMVG(namedtuple("ParamMVG", ["ts", "ti", "tr", "n", "alpha", "Ks"])):
               
    #Recalcul des h0 à partir des param MVG
    def h0(self):
        return -(1/self.alpha)*(((self.ti-self.tr)/(self.ts-self.tr))**(self.n/(1-self.n))-1)**(1/self.n)
    def __repr__(self):
        return "ts{ts}_ti{ti}_tr{tr}_n{n}_alpha{alpha}_Ks{Ks}".format(ts=self.ts,ti=self.ti,tr=self.tr,n=self.n,alpha=self.alpha,Ks=self.Ks)        

# Definition des param gprMax
class ParamGPRMAX :
    pass

def longueur_d_onde(theta, freq, paramMVG, paramGPRMAX):
    """Compute the wavelength for a given water content and frequency (Hz)"""
    eps = CRIM(theta,paramMVG,paramGPRMAX)
    vitesse = 0.3 / math.sqrt(eps) # m/ns
    return vitesse*(10**9)/freq


#%% def des paramètres géométriques
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