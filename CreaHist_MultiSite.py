#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 21:34:02 2020

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
from sklearn.linear_model import LinearRegression
import seaborn as sns


pd.set_option('max_columns', 7)
#%% pathounet
data_path='/home/el/OUT/Codes/Porchet-GPR/OUT_20200722dtrou30_rtrou3_tr5.0/'


X_path='/home/el/Codes/Porchet-GPR/'

hehe=os.getcwd()

if(hehe==X_path):
    print('HEHEHE on est bon')
else:
    os.chdir(X_path)
#%% Reading the folder names
wata='OUT_20200722dtrou30_rtrou3_tr5.0/'
fname=next(os.walk('./'+wata))[1]


#%%

hahat=glob.glob('./Process/*.csv')

#%%
total=pd.DataFrame()
temp=pd.DataFrame()
for filit in hahat:
    temp=pd.read_csv(filit)
    debut=filit.find('_')
    fin=filit.find('_',filit.find('_')+1)
    temp['Lieu']=filit[debut+1:fin]
    total=pd.concat([total, temp])
    #Letter=filit[13]
    
#%%
plt.close('all')
(f1, ax)= plt.subplots(1,1,figsize=(25,15))
sns.violinplot(x="Lieu", y="Ks", data=total,ax=ax,scale="area")
ax.grid()
f1.savefig('./plots/KS_tot.png',format='png')

#sns.swa
    
    