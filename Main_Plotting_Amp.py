#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:58:03 2020

@author: el
"""

import os
import numpy as np
import pandas as pd
from pandas.plotting import parallel_coordinates
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import glob
import io
import numpy as np
from param_acquisition import ParamGPRMAX, ParamMVG, Geometry
from Forward import Forward
from F_extractTWT import F_extractTWT
from F_extract_volumes import F_extract_volumes
import math
import plotly.express as px
from outils import read_parameters, rada_plot


pd.set_option('max_columns', 7)

#%% pathounet
data_path='/home/el/OUT/Codes/Porchet-GPR/OUT_20200722dtrou30_rtrou3_tr5.0/'


X_path='/home/el/Codes/Porchet-GPR/'

hehe=os.getcwd()

if(hehe==X_path):
    print('HEHEHE on est bon')
else:
    os.chdir(X_path)
    
    
    
#%%
foldernama='./OUT_20201121dtrou30_rtrou3_tr5.0/'
fname=next(os.walk(foldernama))[1] 
ouca='Poligny'
#ouca='Bilb'
#ouca='Auffar'
#ouca='Tcherno'
#ouca='Cernay'
fontouney=20
lst=[]

for ii in fname: 

    p=read_parameters(foldernama+ii)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    try:
        temp=F_extractTWT(foldernama+ii)
        Amp=np.genfromtxt(foldernama+ii+'/Amp_EL.csv',delimiter=',',skip_header=1)
        ##
        bibi = 0
        #rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2)/len(TWT_XP))/(max(TWT_XP)-min(TWT_XP))
        #rmsevol=np.sqrt(np.mean((((vol-VOL_XP))**2))/len(VOL_XP))/(max(VOL_XP)-min(VOL_XP))
        ##Si thcerno
        #rmsevol=0
        #rmse=np.sqrt(rmseTwt**2+rmsevol**2)
    except:
        bibi=1
        Amp=np.nan
        #rmseTwt=np.nan
        #rmsevol=np.nan
        #rmse=np.nan
    
    lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,bibi,Amp])
    
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','Converged','Amp']) 

df_params['alpha']=df_params['alpha']*100
df_params['VOL']=df_params['VOL']*0.001

#%%



