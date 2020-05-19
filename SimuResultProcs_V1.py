#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 22:49:07 2020

@author: el
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import io
import numpy as np
from param_acquisition import nT, geometry, ParamMVG, paramGPRMAX, temps

from outils import read_parameters, rada_plot


pd.set_option('max_columns', 7)
#%% pathounet
data_path='/home/el/OUT/Codes/Porchet-GPR/OUTdtrou30_rtrou4_tr5.0/'


X_path='/home/el/Codes/Porchet-GPR/'

hehe=os.getcwd()

if(hehe==X_path):
    print('HEHEHE on est bon')
else:
    os.chdir(X_path)
#%% Reading the folder names
fname=next(os.walk('./OUTdtrou30_rtrou4_tr5.0/'))[1]

#%% pour chaque sous folder, on lit le fichier Params
#s = "Param(a=1, b=2)"



lst=[]

for ii in fname: 

    p=read_parameters('./OUTdtrou30_rtrou4_tr5.0/'+ii)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    try:
        temp = np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', delimiter=',',skip_header=1)
        bibi = 0
        rmse=np.sqrt(np.mean((temp)**2))
    except:
        bibi=1
    
    lst.append([paramMVG.tr,paramMVG.ts,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmse,bibi])
    
df_params=pd.DataFrame(lst,columns=['tr','ts','n','alpha','Ks','RMSE','Converged'])                 


#%% 
plt.close('all')
legendounet=['tr','ts','n','alpha','Ks']
(f1, ax)= plt.subplots(5,5,figsize=(25,15))
ax[0,0].scatter(df_params.tr,df_params.ts,c=df_params.RMSE)
ax[0,0].grid()
ax[0,0].set_xlabel(legendounet[0])
ax[0,0].set_xlabel(legendounet[1])
ax[0,1].scatter(df_params.tr,df_params.n,c=df_params.RMSE)
ax[0,1].grid()
ax[0,1].set_xlabel(legendounet[0])
ax[0,1].set_xlabel(legendounet[2])
ax[1,0].scatter(df_params.tr,df_params.n,c=df_params.RMSE)
ax[1,0].grid()
ax[1,1].scatter(df_params.n,df_params.alpha,c=df_params.RMSE)
ax[1,1].grid()





#%%
pd.plotting.scatter_matrix(df_params)
#%% seqborn
# plt.close('all')
# f1, ax1= plt.subplots(1,1,figsize=(25,15))
# sns.set(style="whitegrid")

#g = sns.PairGrid(df_params, diag_sharey=False)
#g.map_upper(sns.scatterplot)
#g.map_lower(sns.scatterplot)
#g.map_lower(sns.kdeplot, colors="C0")
#g.map_diag(sns.distplot,kde=False)

plt.close('all')
f1, ax1= plt.subplots(1,1,figsize=(25,15))
sns.set(style="ticks")
sns.pairplot(df_params.loc[:,['n','alpha','Ks','ts','tr','Converged']], hue="Converged", markers=["s", "D"])
