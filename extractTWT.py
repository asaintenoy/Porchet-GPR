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


#%% Import du script de définition des param geometrie et GPR max

from param_acquisition import *

#%% import le picking
from picking_radargramme import picking
from outils import read_parameters, rada_plot

#%% parcours le dossier OUTtest
dirName = "OUTtest"
for element in os.listdir(dirName):
    if os.path.isfile(element):
        print("'%s' pas un dossier" % element)
    else:
        ####### Parameters ######
#        filename = 'Parameters'
#        s = io.open(dirName + "/"+ element + "/" + filename).read()    
#        p = eval(s)
        p = read_parameters(dirName + "/" + element + "/")
        print(p)
        print(p[1])
    
        if p[1] == p[1] : #p[1] est ti
            print(p)#Calcul des TWT
            # Definition des paramètres MVG
            paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
            paramMVG.porosity = paramMVG.ts
            
            filename_rada = 'radargram__merged.out'
            
            if (not os.path.isfile(dirName + "/" + element + 'TWT')) :
                print("Pas de TWT dans " + dirName + "/" + element)
                if os.path.isfile(dirName + "/" + element + '/'+ filename_rada) :
                    filename_path =  dirName + "/" + element + '/' + filename_rada
                                
                    print("Picking radargramme dans " + dirName + "/" + element)
                    cas,dt,itmin0,ifenetre,tps_min1,tps_min1_0,tps_min2,tps_min2_0,tps_max,tps_max0,TWT \
                        = picking(filename_path, nT, geometry, paramMVG, paramGPRMAX, temps)  
                    rada_plot(dirName + "/" + element + "/")    
                    # plot des points piqués
                    axe=np.arange(12)
                    plt.plot(axe, (tps_max+tps_max0)/dt,'bo', label='max')        
                    plt.plot(axe, (tps_min1+tps_min1_0)/dt,'ro', label='min1')  
                    plt.plot(axe, (tps_min2+tps_min2_0)/dt,'go', label='min2') 
                    plt.legend(fontsize=20)
                        
                        
                        
