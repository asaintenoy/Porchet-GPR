#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:22:11 2020

@author: el
"""
import matplotlib.pyplot as plt
from param_acquisition import Geometry,ParamMVG, ParamGPRMAX
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import io
import numpy as np
from param_acquisition import ParamGPRMAX, ParamMVG, Geometry
from Forward import Forward
from F_extractTWT import F_extractTWT
from F_extract_volumes import F_extract_volumes

from outils import read_parameters, rada_plot


#%%




TWT_EL=F_extractTWT('./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215_EL/')
VOL_EL=np.genfromtxt('./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215_EL/Volumes_EL.csv',delimiter=',',skip_header=1)


TWT_AS=F_extractTWT('./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215_AS/')
VOL_AS=np.genfromtxt('./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215_AS/Volumes_EL.csv',delimiter=',',skip_header=1)

TWT_S=F_extractTWT('./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215/')
VOL_S=np.genfromtxt('./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215/Volumes_EL.csv',delimiter=',',skip_header=1)



#%%
temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

XX=temps.copy()
XX.insert(0, 0)
plt.close('all')
fig,ax1=plt.subplots(1,1,figsize=(25,15))
ax1.plot(XX,TWT_AS,'.-r',label='TWT_AS(ns)')
ax1.plot(XX,TWT_EL,'--r',label='TWT_EL(ns)')
ax1.plot(XX,TWT_S,'r',label='TWT_S(ns)')
ax1.grid()
ax1.set_xlabel('Minutes')
ax1.set_ylabel('TWT(ns)')
ax2 = ax1.twinx()
ax2.plot(XX,VOL_AS,'.-b',label='Volume_AS')
ax2.plot(XX,VOL_EL,'--b',label='Volume_EL')
ax2.plot(XX,VOL_S,'b',label='Volume_S')
ax2.grid()
ax2.set_ylabel('Volume(mL)')
fig.legend()

#%%
from outils import rada_plot
import h5py
filepath='./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215/'
filename='radargram__merged.out'
f = h5py.File(filepath + filename, 'r')
path = '/rxs/rx1/'
data_S = f['%s%s' % (path, 'Ez')][:,:]

filepath='./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215_EL/'
filename='radargram__merged.out'
f = h5py.File(filepath + filename, 'r')
path = '/rxs/rx1/'
data_EL = f['%s%s' % (path, 'Ez')][:,:]

filepath='./OUTTESRdtrou30_rtrou2_tr10/ts0.35_ti0.07_tr0.011_n5_alpha0.03_Ks0.215_AS/'
filename='radargram__merged.out'
f = h5py.File(filepath + filename, 'r')
path = '/rxs/rx1/'
data_AS = f['%s%s' % (path, 'Ez')][:,:]


(f1, (ax))= plt.subplots(1,3,figsize=(25,15),sharex=True,sharey=True)
# ax[0].xticks(fontsize=20)
# ax[0].yticks(fontsize=20)
ax[0].imshow(data_S[:,:],aspect=0.010)
ax[0].set_xlabel('Homogene',fontsize=20)
ax[1].imshow(data_EL[:,:],aspect=0.010)
ax[1].set_xlabel('EL',fontsize=20)
ax[2].imshow(data_AS[:,:],aspect=0.010)
ax[2].set_xlabel('AS',fontsize=20)
