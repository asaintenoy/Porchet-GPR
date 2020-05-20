#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 14:37:20 2020
@author: el
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


def F_extract_volumes(dirName,temps):
#%% parcours le dossier OUT
    nT=len(temps)
    filename_SWMS = os.path.join(dirName, "SWMS_2D.OUT", "Balance.out")
    if not os.path.isfile(filename_SWMS):
        print("Pas de Balance.out! Probably already simulated")
        
    
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

    return Volume_infiltre
