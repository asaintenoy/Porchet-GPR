#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:48:35 2020
@author: sainteno
"""

#%% Import de librairies
#import itertools
import numpy as np

#import sys
import os
#import h5py
#import math
#import numpy as np
#from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd

#%% Import du script de définition des param geometrie et GPR max
from param_acquisition import nT, geometry, ParamMVG, paramGPRMAX, temps

from outils import read_parameters, rada_plot

#%% parcours le dossier OUT
dir_name = "OUTdtrou30_rtrou4_tr5.0"
PickedVolumes=pd.DataFrame(columns=['Volumes'])

for subdir in os.listdir(dir_name):
    if os.path.isfile(subdir):
        print("'%s' pas un dossier" % subdir)
        continue
    
    print(subdir)
    subdir_name = os.path.join(dir_name, subdir)
    p = read_parameters(subdir_name)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    paramMVG.porosity = paramMVG.ts

    filename_SWMS = os.path.join(subdir_name, "SWMS_2D.OUT", "Balance.out")
    if not os.path.isfile(filename_SWMS):
        print("Pas de Balance.out")
        continue
    
    InFlow=np.zeros(nT+1)
    dVolume=np.zeros(nT)
    Volume_infiltre=np.zeros(nT+1)

    i=0
    fVOL=open(filename_SWMS,"r")

    for ligne in fVOL:
        if 'InFlow' in ligne:
            mots = ligne.split(" ")
            InFlow[i] = float(mots[10])
            i=i+1
    fVOL.close()

    dVolume[0] = temps[0] * InFlow[1]
    Volume_infiltre[1] = Volume_infiltre[0] + dVolume[0]
    
    for i in range(1,nT):
        dVolume[i] = (temps[i] - temps[i-1]) * InFlow[i+1]
        Volume_infiltre[i+1] = Volume_infiltre[i] + dVolume[i]

    PickedVolumes['Volumes']=Volume_infiltre
    PickedVolumes['Volumes'].to_csv(subdir_name + 'Volumes_EL.csv', sep=',', encoding='utf-8', index=None, header=True)
    
#    Pour plotter le volume au cours du temps d'infiltration
#    axe=np.zeros(nT+1)
#    axe[1:nT+1]=temps
#    plt.xlabel("Time (ns)")
#    plt.ylabel("Infiltrated volume of water (cm^3)")
#    plt.plot(axe,Volume_infiltre)        
