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
    
    lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmse,bibi])
    
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSE','Converged'])                 


#%% 
plt.close('all')
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSE','Converged'])
plt.close('all')
legendounet=['ti','ts','n','alpha','Ks']
df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['n']<10.1) ]

(f1, ax)= plt.subplots(5,5,figsize=(25,15))
for ii in range(5):
    for jj in range(5):
        if(ii==jj):
            ax[ii,jj].hist(df_params[legendounet[ii]], weights=np.zeros_like(df_params[legendounet[ii]]) + 1. / size(df_params[legendounet[ii]]))
            ax[ii,jj].set_xlabel(legendounet[ii])
            ax[ii,jj].set_ylabel('Rel Freq.')
            ax[ii,jj].grid()
        else:
            ax[ii,jj].scatter(df_params[legendounet[ii]],df_params[legendounet[jj]],c=df_params.RMSE,cmap = 'jet')
            ax[ii,jj].grid()
            ax[ii,jj].set_xlabel(legendounet[ii])
            ax[ii,jj].set_ylabel(legendounet[jj])
#f1.tight_layout()
        #ax[ii,jj].xaxis.set_label_position('top') 


left  = 0.045  # the left side of the subplots of the figure
right = 0.988    # the right side of the subplots of the figure
bottom = 0.049   # the bottom of the subplots of the figure
top = 0.987      # the top of the subplots of the figure
wspace = 0.224   # the amount of width reserved for blank space between subplots
hspace = 0.290   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

f1.savefig('RMSE.png',format='png')
#plt.close(f1)
#%%
plt.close('all')
legendounet=['tr','ts','n','alpha','Ks']
(f1, ax)= plt.subplots(3,2,figsize=(25,15))
ax[0,0].scatter(df_params['alpha'],df_params.RMSE)
ax[0,0].set_xlabel('alpha')
ax[0,0].grid()
ax[0,1].scatter(df_params['Ks'],df_params.RMSE)
ax[0,1].set_xlabel('Ks')
ax[0,1].grid()
ax[1,1].scatter(df_params['n'],df_params.RMSE)
ax[1,1].set_xlabel('n')
ax[1,1].grid()
ax[2,1].scatter(df_params['ts'],df_params.RMSE)
ax[2,1].set_xlabel('ts')
ax[2,1].grid()
# for ii in range(3):
#     for jj in range(2):
#         ax[ii,jj].scatter(df_params[legendounet[ii]],df_params.RMSE)
#         ax[ii,jj].set_xlabel(legendounet[ii])

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
