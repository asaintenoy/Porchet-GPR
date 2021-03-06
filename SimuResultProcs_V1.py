#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Comparaison pour des données numériques

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
from param_acquisition import ParamGPRMAX, ParamMVG, Geometry
from Forward import Forward
from F_extractTWT import F_extractTWT
from F_extract_volumes import F_extract_volumes

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




#%%

#%% Geometrie
#def des paramètres géométriques
geometry=Geometry()

#Domaine de calcul (en cm)
# largeur
geometry.xmin=0 
geometry.xmax=40
# hauteur (elevation)
geometry.emin=0
geometry.emax = 80
# profondeur du trou en cm
geometry.dtrou = 30
# elevation du fond du trou
geometry.etrou = geometry.emax - geometry.dtrou
 # rayon du trou en cm
geometry.r=4
# hauteur d'eau imposée au fond du trou en cm
geometry.h_eau=5.0
# pas de la maille en cm
geometry.dx = 0.1
#geometry.dx = 1
# profondeur sous le trou (cm) jusqu'où on souhaite un maillage affiné. 
geometry.zaff= 20
#largeur horizontal de la zone affinée (cm)
geometry.waff=20
# elevation de l'affinage
geometry.eaff=geometry.etrou-geometry.zaff 
# contrainte d'angle min pour mesh 
geometry.quality=33
# maximum triangle size  (m*²)
geometry.area=5
# tupple for mesh generation 
geometry.smooth=[1,5]

#%% def de paramètres SWMS
#Temps d'infiltration où à lieu le calcul de chaque trace  (minutes) 
#temps=[1.00, 2.00]
temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

# Temps max de calcul SWMS2D au delà duquel on arrète le calcul (secondes)
tmax_SWMS2D = 600
#tmax_SWMS2D = 10

nT=len(temps)

#%% Definition des param gprMax
paramGPRMAX=ParamGPRMAX()
# Domaine de calcul (cm)
paramGPRMAX.xmin = geometry.xmin
paramGPRMAX.xmax = geometry.xmax
paramGPRMAX.zmin = geometry.emin
paramGPRMAX.zmax = geometry.emax
# Taille des mailles (cm)
paramGPRMAX.dx = 1.0 
# Electrical conductivity of the medium
paramGPRMAX.sigma=0.0000
# Relative dielectric permittivity of water
paramGPRMAX.eps_w=80.1
# Relative dielectric permittivity of PVC
paramGPRMAX.eps_pvc=3
# Relative dielectric permittivity of pure silice
paramGPRMAX.eps_s=2.5
# Ricker signal central frequency (Hz)
paramGPRMAX.wave_freq = 1000e6
# Frequence max du signal EM (Hz)
paramGPRMAX.freq_max = 2.8 * paramGPRMAX.wave_freq
# Distance between hole middle and source (m)
paramGPRMAX.d_emet = 0.18
# Distance between hole middle and receiving antenna (m)
paramGPRMAX.d_recept = 0.22
# param qui raffine le pas spatial (par défaut 10 d'après doc gprmax)
paramGPRMAX.spatial_step = 5
# Trace time window (ns)
paramGPRMAX.time = 30e-9
#time_step_stability_factor (pas utilisé pour le moment...)
paramGPRMAX.fac_dt = 0.2 
#%% pour chaque sous folder, on lit le fichier Params


temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

#Comparé a pour le RMSE
tr = 0.03
# Teneur en eau à saturation
ts = 0.38
# Teneur en eau initiale
ti = 0.09
# Perméabilité à saturation
Ks = 0.2
# param fitting retention n
n = 4
# param fitting retention alpha
alpha = 0.03
pVg=ParamMVG(tr=tr, ts=ts, ti=ti, Ks=Ks, n=n, alpha=alpha)
pVg.porosity = pVg.ts
#[TWT_direct,Vol_direct]=Forward(geometry,pVg,paramGPRMAX,temps,600)
#temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

TWT_direct=F_extractTWT('./OUTdtrou30_rtrou4_tr5.0/ts0.38_ti0.09_tr0.03_n4_alpha0.03_Ks0.2')
Vol_direct=F_extract_volumes('./OUTdtrou30_rtrou4_tr5.0/ts0.38_ti0.09_tr0.03_n4_alpha0.03_Ks0.2',temps)

#%%
lst=[]
for ii in fname: 

    p=read_parameters('./OUTdtrou30_rtrou4_tr5.0/'+ii)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    try:
        #temp = np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', delimiter=',',skip_header=1)
        temp=F_extractTWT('./OUTdtrou30_rtrou4_tr5.0/'+ii)
        vol=np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/Volumes_EL.csv',delimiter=',',skip_header=1)
        bibi = 0
        rmseTwt=np.sqrt(np.mean((temp-TWT_direct)**2))
        rmsevol=np.sqrt(np.mean(((0.001*(vol-Vol_direct))**2)))
        rmse=np.sqrt(rmseTwt**2+rmsevol**2)
    except:
        bibi=1
        rmseTwt=np.nan
        rmsevol=np.nan
        rmse=np.nan
    
    lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmseTwt,rmsevol,rmse,bibi])
    
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSETWT','RMSEVOL','RMSE','Converged'])                 


#%% 
plt.close('all')
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSETWT','RMSEVOL','RMSE','Converged'])
plt.close('all')
legendounet=['ti','ts','n','alpha','Ks']
df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]

(f1, ax)= plt.subplots(5,5,figsize=(25,15))
#cmap = mpl.cm.jet(vmin=0, vmax=1)
#norma = mpl.colors.Normalize(vmin=0, vmax=1)
norm=plt.Normalize(0,2)
for ii in range(5):
    for jj in range(5):
        if(ii==jj):
            ax[ii,jj].hist(df_params[legendounet[ii]], weights=np.zeros_like(df_params[legendounet[ii]]) + 1. / df_params[legendounet[ii]].size)
            ax[ii,jj].set_xlabel(legendounet[ii])
            ax[ii,jj].set_ylabel('Rel Freq.')
            ax[ii,jj].grid()
            
        else:
            sc=ax[ii,jj].scatter(df_params[legendounet[ii]],df_params[legendounet[jj]],c=df_params.RMSE,cmap = 'jet',norm=norm)
            #plt.colorbar(sc,ax=ax[ii,jj])
            ax[ii,jj].grid()
            ax[ii,jj].set_xlabel(legendounet[ii])
            ax[ii,jj].set_ylabel(legendounet[jj])
 
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

f1.savefig('RMSEVOLANDTWT_0306.png',format='png')
#df_params.to_csv('blibalou.csv',sep=',',encoding='utf-8')
#plt.close(f1)
#%%
plt.close('all')
legendounet=['tr','ts','n','alpha','Ks']
(f1, ax)= plt.subplots(2,2,figsize=(25,15))
ax[0,0].scatter(df_params['alpha'],df_params.RMSE)
ax[0,0].set_xlabel('alpha')
ax[0,0].set_ylabel('RMSE')
ax[0,0].grid()
ax[0,1].scatter(df_params['Ks'],df_params.RMSE)
ax[0,1].set_xlabel('Ks')
ax[0,1].set_ylabel('RMSE')
ax[0,1].grid()
ax[1,0].scatter(df_params['n'],df_params.RMSE)
ax[1,0].set_xlabel('n')
ax[1,0].grid()
ax[1,0].set_ylabel('RMSE')

ax[1,1].scatter(df_params['ts'],df_params.RMSE)
ax[1,1].set_xlabel('ts')
ax[1,1].grid()
ax[1,1].set_ylabel('RMSE')
f1.savefig('RMSE-simple.png',format='png')
# for ii in range(3):
#     for jj in range(2):
#         ax[ii,jj].scatter(df_params[legendounet[ii]],df_params.RMSE)
#         ax[ii,jj].set_xlabel(legendounet[ii])


#%%
#df_params.to_csv('3105.csv', sep='\t', encoding='utf-8')
#
#df_params.from_csv('3105.csv')
plt.close('all')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
fig = pyplot.figure()
ax = Axes3D(fig)
gni=ax.scatter(df_params['alpha'],df_params['Ks'],df_params['n'],c=df_params.RMSE,s=abs(df_params.RMSE-df_params.RMSE.max()),cmap = 'jet',vmin=0, vmax=2)
#ax.scatter(df_params['alpha'], df_params['n'], c=df_params.RMSE, zdir='y', zs=1.5)
# ax.plot(y, z, 'g+', zdir='x', zs=-0.5)
# ax.plot(x, y, 'k+', zdir='z', zs=-1.5)

# ax.set_xlim([-0.5, 1.5])
# ax.set_ylim([-0.5, 1.5])
# ax.set_zlim([-1.5, 1.5])



ax.set_xlabel('alpha')
ax.set_ylabel('Ks')
ax.set_zlabel('n')
fig.colorbar(gni)



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
