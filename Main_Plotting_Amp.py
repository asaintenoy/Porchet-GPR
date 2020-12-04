#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:58:03 2020

@author: el
"""

import os
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
    
    
Time_TWT_XP = np.array([0, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 6.00])    
 
#%%
foldernama='./OUT_20200722dtrou30_rtrou3_tr5.0/'
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
        Amp=np.genfromtxt(foldernama+ii+'/AMP_EL.csv',skip_header=1)
        ##
        bibi = 0
        #rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2)/len(TWT_XP))/(max(TWT_XP)-min(TWT_XP))
        #rmsevol=np.sqrt(np.mean((((vol-VOL_XP))**2))/len(VOL_XP))/(max(VOL_XP)-min(VOL_XP))
        ##Si thcerno
        #rmsevol=0
        #rmse=np.sqrt(rmseTwt**2+rmsevol**2)
    except:
        bibi=1
        Amp = np.empty(np.size(Time_TWT_XP))
        Amp[:] = np.nan
        #Amp=np.nan
        #rmseTwt=np.nan
        #rmsevol=np.nan
        #rmse=np.nan
    
    lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,bibi,temp,np.abs(Amp)])
    
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','Converged','TWT','Amp']) 

df_params['alpha']=df_params['alpha']*100
df_params['Meany']=df_params['Amp'].apply(lambda x : np.mean(x))
#df_params['VOL']=df_params['VOL']*0.001

#%% plotting
#cmap = mpl.cm.jet(len(df_params))
plt.close('all')
gni=pd.DataFrame()
gni=df_params.sample(50)
cmap = plt.cm.jet(np.linspace(0,1,len(gni)))
(f2, ax)= plt.subplots(1,2,figsize=(25,15))
for  colo, (index, row) in zip(cmap,gni.iterrows()):
    ax[0].plot(Time_TWT_XP, row['Amp'],c=colo,label=str(np.round(row['Ks'],2))+'-'+str(np.round(row['alpha'],4))+'-'+str(row['n']))
    ax[1].plot(Time_TWT_XP, row['TWT'],c=colo)
    #ax[1,2].plot(Time_TWT_XP, row['VOL'], c=color)
ax[0].grid()
ax[0].legend()
ax[1].grid()
ax[0].set_xlabel('Exp. Time (min)',fontsize=fontouney)
ax[0].set_ylabel('Amp ',fontsize=fontouney)
ax[1].set_xlabel('Exp. Time (min)',fontsize=fontouney)
ax[1].set_ylabel('TWT ',fontsize=fontouney)
#ax[1,2].set_xlabel('Exp. Time (min)',fontsize=fontouney)


#%%En log

plt.close('all')
gni=pd.DataFrame()
gni=df_params.sample(30)
cmap = plt.cm.jet(np.linspace(0,1,len(gni)))
(f2, ax)= plt.subplots(1,1,figsize=(25,15))
for  colo, (index, row) in zip(cmap,gni.iterrows()):
    ax.plot(Time_TWT_XP[2:], np.log(row['Amp'][2:]/np.max(row['Amp'][2:])),c=colo,label=str(np.round(row['Ks'],2))+'-'+str(np.round(row['alpha'],4))+'-'+str(row['n']))
            #ax[1,2].plot(Time_TWT_XP, row['VOL'], c=color)
ax.grid()
ax.legend()

ax.set_xlabel('Exp. Time (min)',fontsize=fontouney)
ax.set_ylabel('Amp ',fontsize=fontouney)


#%% Fitting des droite
from sklearn.linear_model import LinearRegression


plt.close('all')
gni=pd.DataFrame()
gni=df_params.sample(100)
cmap = plt.cm.jet(np.linspace(0,1,len(gni)))
(f2, ax)= plt.subplots(1,2,figsize=(25,15))

for  colo, (index, row) in zip(cmap,gni.iterrows()):
    yy=np.log(row['Amp'][2:]/np.max(row['Amp'][2:]))
    xx=Time_TWT_XP[2:].reshape(-1, 1)
    model = LinearRegression().fit(xx, yy)
    intercept = model.intercept_
    slope = model.coef_

    y_pred = intercept + slope * xx
    ax[0].scatter(xx,yy,60,color=colo,marker='x')
    ax[0].plot(xx,y_pred,color=colo,label=str(np.round(row['Ks'],2))+'-'+str(np.round(row['alpha'],4))+'-'+str(row['n']))
    ax[1].plot(Time_TWT_XP[2:], np.log(row['Amp'][2:]/(np.max(row['Amp'][2:])-np.min(row['Amp'][2:]))),color=colo)
    #np.log(row['Amp'][2:]/np.max(np.abs(row['Amp'][1:])))
            #ax[1,2].plot(Time_TWT_XP, row['VOL'], c=color)
ax[0].grid()
ax[0].legend()
ax[1].grid()
ax[1].legend()
ax[0].set_xlabel('Exp. Time (min)',fontsize=fontouney)
ax[0].set_ylabel('Amp ',fontsize=fontouney)



#%% Appending les params

slopy=[]
intercepty=[]

for  index, row in df_params.iterrows():
    yy=np.log(row['Amp'][2:]/np.max(np.abs(row['Amp'][2:])))
    xx=Time_TWT_XP[2:].reshape(-1, 1)
    if(np.isnan(yy).any()):
        intercept=np.nan
        slope=np.nan
    else:
        model = LinearRegression().fit(xx, yy)
        intercept = model.intercept_
        slope = model.coef_[0]
    
    slopy.append(slope)
    intercepty.append(intercept)
    
#%%
df_params['Amp_slope']=np.array(slopy)
df_params['Amp_intercept']=np.array(intercepty)
df_params['SlopeSign']=np.sign(df_params['Amp_slope'])





#%%
import numpy
import seaborn as sns
plt.close('all')
sns.color_palette("tab10")
sns.pairplot(df_params, hue="species")
#sns.distplot(df_params['Ks'])
g = sns.distplot(data=df_params, x="Ks", hue="SlopeSign", palette="dark", alpha=.6, height=6)




#%%
import seaborn as sns
plt.close('all')
sns.color_palette("tab10")
g = sns.pairplot(df_params[['tr', 'ts', 'ti', 'n', 'alpha', 'Ks','SlopeSign']], hue='SlopeSign',palette='Accent_r')
#g.map_diag(sns.displot)
#g.map_offdiag(sns.scatterplot)
#g.map_upper(sns.scatterplot, s=15)
#g.map_lower(sns.kdeplot)
#g.map_diag(sns.histplot)
#g.add_legend()
#g.savefig("<ultidim.png",dpi=300)




















#%% Histogramme de params
plt.close('all')
legendounet=['ts','n','alpha','Ks']
#df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]

(f1, ax)= plt.subplots(2,2,figsize=(25,15))
#cmap = mpl.cm.jet(vmin=0, vmax=1)
#norma = mpl.colors.Normalize(vmin=0, vmax=1)
norm=plt.Normalize(0,2)
kk=0
for ii in range(2):
    for jj in range(2):
        ax[ii,jj].hist(df_params[df_params['Meany'].isnull()][legendounet[kk]])
        ax[ii,jj].set_xlabel(legendounet[kk])
        ax[ii,jj].set_ylabel('Rel Freq.')
        ax[ii,jj].grid()
        kk=kk+1
            

 
#cbar_ax = f1.add_axes([0.85, 0.15, 0.05, 0.7])
#enfoiros=f1.colorbar(sc, cax=cbar_ax)   
#enfoiros.set_clim(0, 1)
#f1.tight_layout()
        #ax[ii,jj].xaxis.set_label_position('top') 

#plt.colorbar().set_label('Wind speed',rotation=270)
#plt.colorbar(sc)

left  = 0.045  # the left side of the subplots of the figure
right = 0.988    # the right side of the subplots of the figure
bottom = 0.049   # the bottom of the subplots of the figure
top = 0.987      # the top of the subplots of the figure
wspace = 0.224   # the amount of width reserved for blank space between subplots
hspace = 0.290   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
# bordeldenomdedieudemerde=f1.colorbar(sc, ax=ax.ravel().tolist())
# bordeldenomdedieudemerde.set_clim(0, 1)













