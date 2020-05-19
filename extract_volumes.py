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


#%% Import du script de définition des param geometrie et GPR max
from param_acquisition import *

#%% import le picking
#from picking_radargramme import picking
from outils import read_parameters, rada_plot

#%% parcours le dossier OUTtest
dirName = "OUTtest"
for element in os.listdir(dirName):
    if os.path.isfile(element):
        print("'%s' pas un dossier" % element)
        continue
    
    p = read_parameters(dirName + "/" + element + "/")
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    paramMVG.porosity = paramMVG.ts
    
    dirVolume = "SWMS_2D.OUT"
    if not os.path.isfile(dirName + "/" + element + "/" + dirVolume + "/" + "Balance.out"):
        print("Pas de Balance.out")
        continue
    
    InFlow=np.zeros(nT+1)
    dVolume=np.zeros(nT)
    Volume_infiltre=np.zeros(nT+1)

    i=0
    fVOL=open(dirName + "/" + element + "/" + dirVolume + "/" + "Balance.out","r")

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

#    os.chdir(dir + "OUT"+repr(geometry)+"/"+repr(paramMVG))
#    fvolsave=open("Volumes","w")
#    fvolsave.write("""{}\n""".format(Volume_infiltre))
#    fvolsave.close()
#    os.chdir(dir)

    axe=np.zeros(nT+1)
    axe[1:nT+1]=temps
    plt.xlabel("Time (ns)")
    plt.ylabel("Infiltrated volume of water (cm^3)")
    plt.plot(axe,Volume_infiltre)        