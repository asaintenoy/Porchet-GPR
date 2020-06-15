#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Comparaison pour des donnees de terrain

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
import glob
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
wata='OUTdtrou30_rtrou2_tr10.0/'
fname=next(os.walk('./'+wata))[1]

#%% pour config on lit le fichier de data juste pour un
Nama='P1-72-40'
#temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/TWT_'+Nama+'.txt', delimiter=' ')
temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/Auffargis/Twts_Auffar_'+Nama+'.csv', delimiter=',',skip_header=1)
TWT_XP=temp[:,1]
Time_TWT_XP=temp[:,0]
#temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/Volumes_'+Nama+'.txt',delimiter=',')
#VOL_XP_temp=temp[:,1]
#Time_VOL_XP=temp[:,0]
#VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)
temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/Auffargis/Volumes_Auffar_'+Nama+'.csv', delimiter=',',skip_header=1)
VOL_XP=temp[:,1]
Time_VOL_XP=temp[:,0]


#%%
lst=[]
for ii in fname: 

    p=read_parameters('./'+wata+ii)
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    try:
        #temp = np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', delimiter=',',skip_header=1)
        temp=F_extractTWT('./'+wata+ii)
        vol=np.genfromtxt('./'+wata+ii+'/Volumes_EL.csv',delimiter=',',skip_header=1)
        bibi = 0
        rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2))
        rmsevol=np.sqrt(np.mean(((0.001*(vol-VOL_XP))**2)))
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
#df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]

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

f1.savefig('./plots/RMSEVOLANDTWT_'+Nama+'.png',format='png')
#df_params.to_csv('blibalou.csv',sep=',',encoding='utf-8')
#plt.close(f1)
#%% On va regarder la distribution des paramètres pour les 5% des meilleurs modeles

df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False)
df_params_sorted.reset_index(drop=True, inplace=True)
pc=0.1#10percent

df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
df_params_sorted_cut.reset_index(drop=True, inplace=True)

#df_params=df_params[(df_params['RMSE']0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]
plt.close('all')
legendounet=['ts','n','alpha','Ks']
(f2, ax)= plt.subplots(2,2,figsize=(25,15))
kk=0
fontouney=20
for ii in range(2):
    for jj in range(2):
        ax[ii,jj].hist(df_params_sorted_cut[legendounet[kk]], weights=np.zeros_like(df_params_sorted_cut[legendounet[kk]]) + 1. / df_params_sorted_cut[legendounet[kk]].size)
        ax[ii,jj].set_xlabel(legendounet[kk],fontsize=fontouney)
        ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
        ax[ii,jj].grid()        
        kk=kk+1
f2.suptitle(str(100*pc)+'% du meilleur modele - '+Nama, fontsize=fontouney)
f2.savefig('./plots/Histo_'+str(100*pc)+'pc_'+Nama+'.png',format='png')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compa Un par un mais en boucle
#hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/twt*.txt')
#hahav=glob.glob('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/vol*.txt')

# hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/TWT*.txt')
# hahav=glob.glob('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/Vol*.txt')
hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/Auffargis/Twts_*.csv')
hahav=glob.glob('/home/el/Data/Compil_data-Kriterres/Auffargis/Volumes_*.csv')


lst=[]
#%% Reading the folder names
foldernama='./OUTdtrou30_rtrou2_tr10.0/'
fname=next(os.walk(foldernama))[1]

#%%
#ouca='Poligny'
#ouca='Bilbo'
ouca='Auffargis_30'
fontouney=20
for filit,filiv in zip(hahat,hahav):
    lst=[]
    #Nama=filit[-11:-4] #Poligny
    #Nama=filit[-12:-4]#Bilbo
    Nama=filit[-12:-4]#Auffargis
    #temp=np.genfromtxt(filit, delimiter=' ')
    temp=np.genfromtxt(filit, delimiter=',',skip_header=1)
    TWT_XP=temp[:,1]
    Time_TWT_XP=temp[:,0]
    temp=np.genfromtxt(filiv,delimiter=',',skip_header=1)
    VOL_XP_temp=temp[:,1]
    Time_VOL_XP=temp[:,0]
    VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)
    
    for ii in fname: 
    
        p=read_parameters(foldernama+ii)
        paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
        try:
            #temp = np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', delimiter=',',skip_header=1)
            temp=F_extractTWT(foldernama+ii)
            vol=np.genfromtxt(foldernama+ii+'/Volumes_EL.csv',delimiter=',',skip_header=1)
            bibi = 0
            rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2))
            rmsevol=np.sqrt(np.mean(((0.001*(vol-VOL_XP))**2)))
            rmse=np.sqrt(rmseTwt**2+rmsevol**2)
        except:
            bibi=1
            rmseTwt=np.nan
            rmsevol=np.nan
            rmse=np.nan
        
        lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmseTwt,rmsevol,rmse,bibi,(temp),(vol)])
        
    df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSETWT','RMSEVOL','RMSE','Converged','TWT','VOL']) 
    
    
    
################Sensibility plot
    # plt.close('all')
    # legendounet=['ti','ts','n','alpha','Ks']
    # #df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]
    
    # (f1, ax)= plt.subplots(5,5,figsize=(25,15))
    # #cmap = mpl.cm.jet(vmin=0, vmax=1)
    # #norma = mpl.colors.Normalize(vmin=0, vmax=1)
    # norm=plt.Normalize(0,2)
    # for ii in range(5):
    #     for jj in range(5):
    #         if(ii==jj):
    #             ax[ii,jj].hist(df_params[legendounet[ii]], weights=np.zeros_like(df_params[legendounet[ii]]) + 1. / df_params[legendounet[ii]].size)
    #             ax[ii,jj].set_xlabel(legendounet[ii],fontsize=fontouney)
    #             ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
    #             ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
    #             ax[ii,jj].grid() 
                
    #         else:
    #             sc=ax[ii,jj].scatter(df_params[legendounet[ii]],df_params[legendounet[jj]],c=df_params.RMSE,cmap = 'jet',norm=norm)
    #             #plt.colorbar(sc,ax=ax[ii,jj])
    #             ax[ii,jj].grid()
    #             ax[ii,jj].set_xlabel(legendounet[ii],fontsize=fontouney)
    #             ax[ii,jj].set_ylabel(legendounet[jj],fontsize=fontouney)
    #             ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)

     
    
    # left  = 0.045  # the left side of the subplots of the figure
    # right = 0.988    # the right side of the subplots of the figure
    # bottom = 0.049   # the bottom of the subplots of the figure
    # top = 0.987      # the top of the subplots of the figure
    # wspace = 0.224   # the amount of width reserved for blank space between subplots
    # hspace = 0.290   # the amount of height reserved for white space between subplots
    # plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    # f1.savefig('./plots/RMSEVOLANDTWT_'+ouca+'-'+Nama+'.png',format='png')
    
####################    
################### histogram
    df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False)
    df_params_sorted.reset_index(drop=True, inplace=True)
    pc=0.1#10percent
    
    df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
    df_params_sorted_cut.reset_index(drop=True, inplace=True)

    #df_params=df_params[(df_params['RMSE']0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]
    plt.close('all')
    legendounet=['ts','n','alpha','Ks']
    (f2, ax)= plt.subplots(2,3,figsize=(25,15))
    kk=0
    
    for ii in range(2):
        for jj in range(2):
            ax[ii,jj].hist(df_params_sorted_cut[legendounet[kk]], weights=np.zeros_like(df_params_sorted_cut[legendounet[kk]]) + 1. / df_params_sorted_cut[legendounet[kk]].size)
            ax[ii,jj].set_xlabel(legendounet[kk],fontsize=fontouney)
            ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
            ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
            ax[ii,jj].grid() 
            ax[ii,jj].set_ylim([0, 0.3])
            kk=kk+1
    
    cmap = mpl.cm.jet((df_params_sorted_cut['RMSE'].values - df_params_sorted_cut['RMSE'].min())/\
    (df_params_sorted_cut['RMSE'].max() - df_params_sorted_cut['RMSE'].min()))    
    for  color, (index, row) in zip(cmap,df_params_sorted_cut.iterrows()):
        ax[0,2].plot(Time_TWT_XP, row['TWT'], c=color)
        ax[1,2].plot(Time_TWT_XP, row['VOL'], c=color)
        ax[0,2].grid()
        ax[1,2].grid()
    ax[0,2].plot(Time_TWT_XP, TWT_XP, marker = '*', markersize = 10, mfc = 'k')
    ax[1,2].plot(Time_TWT_XP, VOL_XP, marker = '*', markersize = 10, mfc = 'k')    
    f2.suptitle(str(100*pc)+'% du meilleur modele '+ ouca+'-'+Nama, fontsize=fontouney)
    f2.savefig('./plots/Histo_'+str(100*pc)+'pc_'+ouca+'-'+Nama+'.png',format='png')
###################















#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compa total les sommant tous

#%% on lit les fichiers de data
#filiname=next(os.walk('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/*.txt'))[2]



hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/TWT*.txt')
hahav=glob.glob('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/Vol*.txt')
lst=[]
#%% Reading the folder names
fname=next(os.walk('./OUTdtrou30_rtrou4_tr5.0/'))[1]
#%%
for filit,filiv in zip(hahat,hahav):
    #Nama=fili[-16:-4]
    temp=np.genfromtxt(filit, delimiter=' ')
    TWT_XP=temp[:,1]
    Time_TWT_XP=temp[:,0]
    temp=np.genfromtxt(filiv,delimiter=',')
    VOL_XP_temp=temp[:,1]
    Time_VOL_XP=temp[:,0]
    VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)
    
    for ii in fname: 
    
        p=read_parameters('./OUTdtrou30_rtrou4_tr5.0/'+ii)
        paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
        try:
            #temp = np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', delimiter=',',skip_header=1)
            temp=F_extractTWT('./OUTdtrou30_rtrou4_tr5.0/'+ii)
            vol=np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/Volumes_EL.csv',delimiter=',',skip_header=1)
            bibi = 0
            rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2))
            rmsevol=np.sqrt(np.mean(((0.001*(vol-VOL_XP))**2)))
            rmse=np.sqrt(rmseTwt**2+rmsevol**2)
        except:
            bibi=1
            rmseTwt=np.nan
            rmsevol=np.nan
            rmse=np.nan
        
        lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmseTwt,rmsevol,rmse,bibi])
        
    df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSETWT','RMSEVOL','RMSE','Converged'])  

#%%
df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False)
df_params_sorted.reset_index(drop=True, inplace=True)
pc=0.1#10percent

df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
df_params_sorted_cut.reset_index(drop=True, inplace=True)

#df_params=df_params[(df_params['RMSE']0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]
plt.close('all')
legendounet=['ts','n','alpha','Ks']
(f2, ax)= plt.subplots(2,2,figsize=(25,15))
kk=0
fontouney=20
for ii in range(2):
    for jj in range(2):
        ax[ii,jj].hist(df_params_sorted_cut[legendounet[kk]], weights=np.zeros_like(df_params_sorted_cut[legendounet[kk]]) + 1. / df_params_sorted_cut[legendounet[kk]].size)
        ax[ii,jj].set_xlabel(legendounet[kk],fontsize=fontouney)
        ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
        ax[ii,jj].grid()        
        kk=kk+1
f2.suptitle(str(100*pc)+'% du meilleur modele - Total-Poligny', fontsize=fontouney)
f2.savefig('./plots/Histo_'+str(100*pc)+'pc_Total-Poligny.png',format='png')










