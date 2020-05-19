#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:19:41 2020
@author: sainteno
"""

#%% Import de librairies
#import itertools
import numpy as np

#import sys
import os
import h5py
#import math
#import numpy as np
#from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd


#%% Import du script de définition des param geometrie et GPR max

from param_acquisition import nT, geometry, ParamMVG, paramGPRMAX, temps

#%% import le picking
from picking_radargramme import picking
from outils import read_parameters, rada_plot

#%% parcours le dossier OUT
dir_name = "OUTdtrou30_rtrou4_tr5.0"
PickedTWT=pd.DataFrame(columns=['TWT(ns)'])
for subdir in os.listdir(dir_name):
    if os.path.isfile(subdir):
        print("'%s' pas un dossier" % subdir)
        continue
    
    print(subdir)
    subdir_name = os.path.join(dir_name, subdir)
    p = read_parameters(subdir_name)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    paramMVG.porosity = paramMVG.ts

    filename_TWT = os.path.join(subdir_name, 'TWT_EL.csv')
    if os.path.isfile(filename_TWT) :
        continue
 
    filename_radar =  os.path.join(subdir_name, 'radargram__merged.out')
    
    if not os.path.isfile(filename_radar) :
        print("pas de radargramme dans " + subdir_name)
        continue

    print("Picking radargramme dans " + subdir_name)
    
    cas,dt,itmin0,ifenetre,tps_min1,tps_min1_0,tps_min2,tps_min2_0,tps_max,tps_max0,TWT \
        = picking(filename_radar, nT, geometry, paramMVG, paramGPRMAX, temps)  
        
    
    PickedTWT['TWT(ns)']=TWT
    PickedTWT['TWT(ns)'].to_csv(filename_TWT, sep=',', encoding='utf-8', index=None, header=True)

    #rada_plot(subdir_name)    
    # plot des points piqués
#    axe=np.arange(12)
#    plt.plot(axe, (tps_max+tps_max0)/dt,'bo', label='max')        
#    plt.plot(axe, (tps_min1+tps_min1_0)/dt,'ro', label='min1')  
#    plt.plot(axe, (tps_min2+tps_min2_0)/dt,'go', label='min2') 
#    plt.legend(fontsize=20)


