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
def F_extractTWT(pathounet):
    temp = np.genfromtxt(pathounet+'/TWT_EL.csv', delimiter=',',skip_header=1)
    #try:
    #    temp = np.genfromtxt(pathounet+'/TWT_EL.csv', delimiter=',',skip_header=1)
    # except:
    #    temp = np.empty((11))
    #    temp[:]=np.nan
       
    return temp
       