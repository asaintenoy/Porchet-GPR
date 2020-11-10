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
#%% Reading the folder names
wata='OUT_20200722dtrou30_rtrou3_tr5.0/'
fname=next(os.walk('./'+wata))[1]

#%% pour config on lit le fichier de data juste pour un
Nama='BilbP1_0_30'
#temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/TWT_'+Nama+'.txt', delimiter=' ')
#temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/Auffargis/Twts_Auffar_'+Nama+'.csv', delimiter=',',skip_header=1)

#hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/twt*.txt')
#hahav=glob.glob('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/vol*.txt')

temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/twts_'+Nama+'.csv', delimiter=',',skip_header=0)
TWT_XP=temp[:,1]
Time_TWT_XP=temp[:,0]
#temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/Volumes_'+Nama+'.txt',delimiter=',')
#VOL_XP_temp=temp[:,1]
#Time_VOL_XP=temp[:,0]
#VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)volumes_BilbP1_0_30.txt
temp=np.genfromtxt('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/volumes_'+Nama+'.csv', delimiter=',',skip_header=0)
VOL_XP_temp=temp[:,1]
Time_VOL_XP=temp[:,0]
VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)

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
legendounet=['ti','ts','n','alpha','Ks']
gnifalou=df_params['ti'].unique()
shortl=['tr','ts','n','alpha','Ks']



fontouney=20
for kk in gnifalou:
    (f1, ax)= plt.subplots(5,5,figsize=(25,15))
    norm=plt.Normalize(0,2)
    df_params_cutted=df_params[(df_params['ti']==round(float(kk),2)) & (df_params['n']>=5)].copy()
    for ii in range(5):
        for jj in range(5):
            if(ii==jj):
                ax[ii,jj].hist(df_params_cutted[shortl[ii]], weights=np.zeros_like(df_params_cutted[shortl[ii]]) + 1. / df_params_cutted[shortl[ii]].size)
                ax[ii,jj].set_xlabel(shortl[ii])
                ax[ii,jj].set_ylabel('Rel Freq.')
                ax[ii,jj].grid()
                
            else:
                ax[ii,jj].scatter(df_params_cutted[shortl[ii]],df_params_cutted[shortl[jj]],c=df_params_cutted.RMSE,cmap = 'jet',norm=norm)
                #plt.colorbar(sc,ax=ax[ii,jj])
                ax[ii,jj].grid()
                ax[ii,jj].set_xlabel(shortl[ii])
                ax[ii,jj].set_ylabel(shortl[jj])

    f1.suptitle(' ti = '+str(kk)+'n>=5', fontsize=fontouney)
    f1.savefig('./plots/RMSETWT_Test_ti='+str(kk)+'n>=5.png',format='png')


# for ii in range(2):
#     for jj in range(2):
#             ax[ii,jj].scatter(df_params_cutted[shortl[ii]],df_params_cutted[shortl[jj]],c=df_params_cutted.RMSE,cmap = 'jet',norm=norm)
#             #plt.colorbar(sc,ax=ax[ii,jj])
#             ax[ii,jj].grid()
#             ax[ii,jj].set_xlabel(shortl[ii])
#             ax[ii,jj].set_ylabel(shortl[jj])
#             cc=cc+1
    

#%% 
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

f1.savefig('./plots/RMSETWT_'+Nama+'.png',format='png')
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

#Poligny
#hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/190527-Poligny/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/TWT*.csv')

#CDC
hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/twt*.csv')

#Auffargis
#hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/Auffargis/Twts_Auffar*.csv')

#Tcherno
#hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/Tcherno/*-20-30.csv')

#CErnay 2020
#hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/Cernay/Twts_Cernay_*.csv')




lst=[]
#%% Reading the folder names
#bilbo+Poligny
#foldernama='./OUTdtrou35_rtrou4_tr5.0/'
#foldernama='./OUTdtrou30_rtrou4_tr5.0/'

#Auffar
#foldernama='./OUTdtrou30_rtrou2_tr10.0/'
#tcherno
#foldernama='./OUT_chernodtrou30_rtrou2_5_tr10/'

#CErnay
#polignyet bilbo nouveau
foldernama='./OUT_20200722dtrou30_rtrou3_tr5.0/'

fname=next(os.walk(foldernama))[1]





#%%
#ouca='Poligny'
ouca='Bilb'
#ouca='Auffar'
#ouca='Tcherno'
#ouca='Cernay'
fontouney=20
for filit in hahat:
    lst=[]
    fin=filit.find('.csv',-4)   
    #guili=filit.find('/',-25)#Auffargis
    #guili=filit.find('/',-20)#CErnay
    guili=filit.find('/',-40)#Bilbo
    #guili=filit.find('/',-40)#Poligny
    debut=filit.find(ouca,guili)+len(ouca)
    
    #Nama=filit[debut:fin]
    #Nama=filit[-11:-4] #Poligny
    Nama=filit[-12:-4]#Bilbo
    #Nama=filit[-12:-4]#Auffargis
    #Nama=filit[-7:-4]#CErnay
    #temp=np.genfromtxt(filit, delimiter=',')
    temp=np.genfromtxt(filit, delimiter=',',skip_header=1)
    #pour poiligny et autre merde a la con resempler
    TWT_XP_RS = np.array([0, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 6.00])    
    #TWT_XP_1=temp[:,1]
    #Time_TWT_XP_1=temp[:,0]
    TWT_XP=np.interp(TWT_XP_RS,temp[:,0],temp[:,1])
    Time_TWT_XP=TWT_XP_RS



    #temp=np.genfromtxt(filit[0:guili]+'/'+'volumes_'+ouca+Nama+'.csv',delimiter=',',skip_header=1)            
    temp=np.genfromtxt(filit[0:guili]+'/'+'Volumes_'+ouca+'_'+Nama+'.csv',delimiter=',',skip_header=1)  
    VOL_XP_temp=temp[:,1]
    Time_VOL_XP=temp[:,0]
    VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)
    
    #Si tcherno
    

    
    for ii in fname: 
    
        p=read_parameters(foldernama+ii)
        paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
        try:
            temp=F_extractTWT(foldernama+ii)
            vol=np.genfromtxt(foldernama+ii+'/Volumes_EL.csv',delimiter=',',skip_header=1)
            ##
            bibi = 0
            rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2)/len(TWT_XP))/(max(TWT_XP)-min(TWT_XP))
            rmsevol=np.sqrt(np.mean((((vol-VOL_XP))**2))/len(VOL_XP))/(max(VOL_XP)-min(VOL_XP))
            ##Si thcerno
            #rmsevol=0
            rmse=np.sqrt(rmseTwt**2+rmsevol**2)
        except:
            bibi=1
            rmseTwt=np.nan
            rmsevol=np.nan
            rmse=np.nan
        
        lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmseTwt,rmsevol,rmse,bibi,(temp),(vol)])
        
    df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSETWT','RMSEVOL','RMSE','Converged','TWT','VOL']) 
    
    df_params['alpha']=df_params['alpha']*100
    df_params['VOL']=df_params['VOL']*0.001
    
    
    
    pc=0.01#percent
    
################Sensibility plot

    df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False,ascending=True)
    df_params_sorted.reset_index(drop=True, inplace=True)
    
    
    df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
    df_params_sorted_cut.reset_index(drop=True, inplace=True)
    


    plt.close('all')
    legendounet=['ti','ts','n','alpha','Ks']
    #df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]
    
    (f1, ax)= plt.subplots(5,5,figsize=(25,15))
    #cmap = mpl.cm.jet(vmin=0, vmax=1)
    #norma = mpl.colors.Normalize(vmin=0, vmax=1)
    #maxi=df_params_sorted_cut['RMSE'][len(df_params_sorted_cut)-1]
    #mini=df_params_sorted_cut['RMSE'][0]
    #maxi=np.log10(1.5)
    #mini=np.log10(10**(-1))
    maxi=0.03#0.06 pour 10% 0.03 pour 1%
    mini=0.01
    #df_params['RMSE'] = df_params['RMSE'].apply(np.log10)
    norm=plt.Normalize(mini,maxi)
    df_params_sorted_cut.sort_values(by=['RMSE'],inplace=True,ascending=False)
    for ii in range(5):
        for jj in range(5):
            if(ii==jj):
                ax[ii,jj].hist(df_params_sorted_cut[legendounet[ii]], weights=np.zeros_like(df_params_sorted_cut[legendounet[ii]]) + 1. / df_params_sorted_cut[legendounet[ii]].size)
                ax[ii,jj].set_xlabel(legendounet[ii],fontsize=fontouney)
                ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
                ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
                ax[ii,jj].grid() 
                
            else:
                sc=ax[ii,jj].scatter(df_params_sorted_cut[legendounet[ii]],df_params_sorted_cut[legendounet[jj]],c=df_params_sorted_cut.RMSE,cmap = 'jet',norm=norm)
                #sc=ax[ii,jj].scatter(df_params_sorted_cut[legendounet[ii]],df_params_sorted_cut[legendounet[jj]],c=df_params_sorted_cut.RMSE,cmap = 'jet',norm=norm)
                #plt.colorbar(sc,ax=ax[ii,jj])
                ax[ii,jj].grid()
                ax[ii,jj].set_xlabel(legendounet[ii],fontsize=fontouney)
                ax[ii,jj].set_ylabel(legendounet[jj],fontsize=fontouney)
                ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)

     
    
    left  = 0.045  # the left side of the subplots of the figure
    right = 0.988    # the right side of the subplots of the figure
    bottom = 0.049   # the bottom of the subplots of the figure
    top = 0.987      # the top of the subplots of the figure
    wspace = 0.224   # the amount of width reserved for blank space between subplots
    hspace = 0.290   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    f1.savefig('./plots/'+ouca+'/NormalizedRMSEVOLANDTWT_'+str(np.round(100*pc))+'pc_'+ouca+Nama+'.png',format='png')
    
####################   Parallel coordinates

    df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False,ascending=True)
    df_params_sorted.reset_index(drop=True, inplace=True)
    
    df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
    df_params_sorted_cut.reset_index(drop=True, inplace=True)
    cmap = mpl.cm.jet((df_params_sorted_cut['RMSE'].values - df_params_sorted_cut['RMSE'].min())/\
    (df_params_sorted_cut['RMSE'].max() - df_params_sorted_cut['RMSE'].min()))                     
    # parallel_coordinates(df_params_sorted_cut,class_column='RMSE',color=cmap,cols=['ti','ts','alpha', 'n', 'Ks'])


    fig=px.parallel_coordinates(df_params_sorted_cut[['ti','ts','alpha', 'n', 'Ks','RMSE']], color='RMSE', labels={'ti','ts','alpha', 'n', 'Ks'},\
                            color_continuous_scale=px.colors.diverging.Tealrose)
    #fig.show()
    fig.write_html('./plots/'+ouca+'/NormalizedRMSEParallelC_'+str(np.round(100*pc))+'pc_'+ouca+Nama+'.html')
 
################### histogram
    df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False,ascending=True)
    df_params_sorted.reset_index(drop=True, inplace=True)

    
    df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
    df_params_sorted_cut.reset_index(drop=True, inplace=True)
    
    #df_params_sorted_cut=df_params_sorted[df_params_sorted['RMSE']<0.1]
    #df_params_sorted_cut.reset_index(drop=True, inplace=True)

    plt.close('all')
    legendounet=['ts','n','alpha','Ks']
    (f2, ax)= plt.subplots(2,3,figsize=(25,15))
    kk=0
    df_params_sorted_cut.sort_values(by=['RMSE'],inplace=True,ascending=False)

    for ii in range(2):
        for jj in range(2):
            ax[ii,jj].hist(df_params_sorted_cut[legendounet[kk]], weights=np.zeros_like(df_params_sorted_cut[legendounet[kk]]) + 1. / df_params_sorted_cut[legendounet[kk]].size)
            ax[ii,jj].set_xlabel(legendounet[kk]+' Optim= '+str(df_params_sorted_cut[legendounet[kk]][0]),fontsize=fontouney)
            ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
            #ax[ii,jj].plot(df_params_sorted_cut[legendounet[kk]][0], 0.15, marker = '*', markersize = 10, mfc = 'k')
            ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
            ax[ii,jj].grid() 
            ax[ii,jj].set_ylim([0, 0.3])
            kk=kk+1




    # #tcherno    
    # ii=1
    # jj=2
    # legendounet=['ti']
    # kk=0
    # ax[ii,jj].hist(df_params_sorted_cut[legendounet[kk]], weights=np.zeros_like(df_params_sorted_cut[legendounet[kk]]) + 1. / df_params_sorted_cut[legendounet[kk]].size)
    # ax[ii,jj].set_xlabel(legendounet[kk]+' Optim= '+str(df_params_sorted_cut[legendounet[kk]][0]),fontsize=fontouney)
    # ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
    # #ax[ii,jj].plot(df_params_sorted_cut[legendounet[kk]][0], 0.15, marker = '*', markersize = 10, mfc = 'k')
    # ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
    # ax[ii,jj].grid() 
    # ax[ii,jj].set_ylim([0, 0.3])
    
    cmap = mpl.cm.jet((df_params_sorted_cut['RMSE'].values - df_params_sorted_cut['RMSE'].min())/\
    (df_params_sorted_cut['RMSE'].max() - df_params_sorted_cut['RMSE'].min()))    
    for  color, (index, row) in zip(cmap,df_params_sorted_cut.iterrows()):
        ax[0,2].plot(Time_TWT_XP, row['TWT'], c=color)
        ax[1,2].plot(Time_TWT_XP, row['VOL'], c=color)
    ax[0,2].grid()
    ax[1,2].grid()
    ax[0,2].plot(Time_TWT_XP, TWT_XP, marker = '+', markersize = 10, mfc = 'k')
    ax[0,2].plot(Time_TWT_XP,df_params_sorted_cut['TWT'][0], marker = '*', markersize = 10, mfc = 'k')
    ax[1,2].plot(Time_TWT_XP, 0.001*VOL_XP, marker = '+', markersize = 10, mfc = 'k') 
    ax[1,2].plot(Time_TWT_XP,df_params_sorted_cut['VOL'][0], marker = '*', markersize = 10, mfc = 'k')
    ax[0,2].tick_params(axis='both', which='major', labelsize=fontouney)
    ax[1,2].tick_params(axis='both', which='major', labelsize=fontouney)
    ax[0,2].set_ylabel('TWT(ns)',fontsize=fontouney)
    ax[0,2].set_xlabel('Exp. Time (min)',fontsize=fontouney)
    ax[1,2].set_ylabel('Vol. (l)',fontsize=fontouney)
    ax[1,2].set_xlabel('Exp. Time (min)',fontsize=fontouney)
    
    
    f2.suptitle(str(np.round(100*pc))+'% des modèles '+ ouca+'-'+Nama+'\n Nombre de modeles: '+str(len(df_params_sorted_cut))+', RMSE: ['+str(np.round(df_params_sorted_cut['RMSE'][0],4))+';'+str(np.round(df_params_sorted_cut['RMSE'][len(df_params_sorted_cut)-1],4))+'] (ns)', fontsize=fontouney)
    f2.savefig('./plots/'+ouca+'/Histo_'+'pc_'+str(round(100*pc))+'_'+ouca+'-'+Nama+'.png',format='png')
##################

#############Sensibility plot mais en coupant Ti

    # plt.close('all')
    # legendounet=['ti','ts','n','alpha','Ks']
    # gnifalou=df_params['ti'].unique()
    # shortl=['tr','ts','n','alpha','Ks']
    
    # fontouney=20
    # for kk in gnifalou:
    #     (f1, ax)= plt.subplots(5,5,figsize=(25,15))
    #     norm=plt.Normalize(0,2)
    #     df_params_cutted=df_params[df_params['ti']==round(float(kk),2)].copy()
    #     for ii in range(5):
    #         for jj in range(5):
    #             if(ii==jj and df_params_cutted.empty==False):
    #                 ax[ii,jj].hist(df_params_cutted[shortl[ii]], weights=np.zeros_like(df_params_cutted[shortl[ii]]) + 1. / df_params_cutted[shortl[ii]].size)
    #                 ax[ii,jj].set_xlabel(shortl[ii])
    #                 ax[ii,jj].set_ylabel('Rel Freq.')
    #                 ax[ii,jj].grid()
                    
    #             else:
    #                 ax[ii,jj].scatter(df_params_cutted[shortl[ii]],df_params_cutted[shortl[jj]],c=df_params_cutted.RMSE,cmap = 'jet',norm=norm)
    #                 #plt.colorbar(sc,ax=ax[ii,jj])
    #                 ax[ii,jj].grid()
    #                 ax[ii,jj].set_xlabel(shortl[ii])
    #                 ax[ii,jj].set_ylabel(shortl[jj])   
    #     f1.suptitle(' ti = '+str(kk), fontsize=fontouney)
    #     f1.savefig('./plots/'+ouca+'/RMSEVOLANDTWT_'+ouca+Nama+'_ti='+str(kk)+'.png',format='png')



##### RMSE en fonction de param

    # plt.close('all')
    # legendounet=['ts','n','alpha','Ks']
    # (f2, ax)= plt.subplots(2,2,figsize=(25,15))
    # kk=0
    # norm=plt.Normalize(0,2)
    # for ii in range(2):
    #     for jj in range(2):
    #         ax[ii,jj].scatter(df_params[legendounet[kk]],df_params['RMSE'],c=df_params.RMSE,cmap = 'jet',norm=norm)
    #         ax[ii,jj].set_xlabel(legendounet[kk]+' Optim= '+str(df_params_sorted_cut[legendounet[kk]][0]),fontsize=fontouney)
    #         ax[ii,jj].set_ylabel('RMSE',fontsize=fontouney)
    #         #ax[ii,jj].plot(df_params_sorted_cut[legendounet[kk]][0], 0.15, marker = '*', markersize = 10, mfc = 'k')
    #         ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
    #         ax[ii,jj].grid() 
    #         ax[ii,jj].set_ylim([0, 5])
    #         kk=kk+1
    # f2.savefig('./plots/'+ouca+'/RMSE_Param'+ouca+'-'+Nama+'.png',format='png')




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





# #%% Echelle de couleur
# ouca='Tcherno'
# fontouney=20
# for filit in hahat:
#     lst=[]
#     fin=filit.find('.csv',-4)   
#     guili=filit.find('/',-40)
#     debut=filit.find(ouca,guili)+len(ouca)
    
#     Nama=filit[debut:fin]
#     #Nama=filit[-11:-4] #Poligny
#     #Nama=filit[-12:-4]#Bilbo
#     #Nama=filit[-12:-4]#Auffargis
#     #temp=np.genfromtxt(filit, delimiter=',')
#     temp=np.genfromtxt(filit, delimiter=',',skip_header=1)
#     TWT_XP=temp[:,1]
#     Time_TWT_XP=temp[:,0]



#     # temp=np.genfromtxt(filit[0:guili]+'/'+'Volumes_'+ouca+Nama+'.csv',delimiter=',',skip_header=1)            
#     # VOL_XP_temp=temp[:,1]
#     # Time_VOL_XP=temp[:,0]
#     # VOL_XP=np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)
    
#     for ii in fname: 
    
#         p=read_parameters(foldernama+ii)
#         paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
#         try:
#             #temp = np.genfromtxt('./OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', delimiter=',',skip_header=1)
#             temp=F_extractTWT(foldernama+ii)
#             #vol=np.genfromtxt(foldernama+ii+'/Volumes_EL.csv',delimiter=',',skip_header=1)
#             vol=np.nan
#             bibi = 0
#             rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2))
#             #rmsevol=np.sqrt(np.mean(((0.001*(vol-VOL_XP))**2)))
#             rmsevol=0
#             rmse=np.sqrt(rmseTwt**2+rmsevol**2)
#         except:
#             bibi=1
#             rmseTwt=np.nan
#             rmsevol=np.nan
#             rmse=np.nan
        
#         lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,rmseTwt,rmsevol,rmse,bibi,(temp),(vol)])
        
#     df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','RMSETWT','RMSEVOL','RMSE','Converged','TWT','VOL']) 
    
#     df_params['alpha']=df_params['alpha']*100
#     df_params['VOL']=df_params['VOL']*0.001
# ################Sensibility plot

#     df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False,ascending=False)
#     df_params_sorted.reset_index(drop=True, inplace=True)
#     pc=1#percent
    
#     df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
#     df_params_sorted_cut.reset_index(drop=True, inplace=True)
    
#     # df_params_sorted_cut=df_params_sorted[df_params_sorted['RMSE']<0.1]
#     # df_params_sorted_cut.reset_index(drop=True, inplace=True)
#     df_params=pd.DataFrame()
#     df_params=df_params_sorted_cut.copy()

#     plt.close('all')
#     legendounet=['ti','ts','n','alpha','Ks']
#     #df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]
    
#     (f1, ax)= plt.subplots(1,1,figsize=(25,15))
#     #cmap = mpl.cm.jet(vmin=0, vmax=1)
#     #norma = mpl.colors.Normalize(vmin=0, vmax=1)
#     maxi=df_params_sorted_cut['RMSE'][len(df_params)-1]
#     mini=df_params_sorted_cut['RMSE'][0]
#     maxi=np.log10(1)
#     mini=np.log10(0.05)
#     #maxi=np.log10(2)#0.6
#     #mini=np.log10(0.1)
#     #df_params['RMSE'] = df_params['RMSE'].apply(np.log10)
#     norm=plt.Normalize(mini,maxi)
#     ii=3
#     jj=4
#     sc=ax.scatter(df_params_sorted_cut[legendounet[ii]],df_params_sorted_cut[legendounet[jj]],c=df_params_sorted_cut.RMSE,cmap = 'jet',norm=norm)
#     cbar=plt.colorbar(sc)
#     cbar.ax.tick_params(labelsize=fontouney)

    
#     #f1.savefig('./plots/'+ouca+'/RMSEVOLANDTWT_'+str(np.round(100*pc))+'pc_'+ouca+Nama+'.png',format='png')




