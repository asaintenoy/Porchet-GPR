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
import math

from outils import read_parameters, rada_plot


pd.set_option('max_columns', 7)
#%% pathounet
data_path='/home/el/OUT/Codes/Porchet-GPR/OUTdtrou35_rtrou4_tr5.0/'


X_path='/home/el/Codes/Porchet-GPR/'

hehe=os.getcwd()

if(hehe==X_path):
    print('HEHEHE on est bon')
else:
    os.chdir(X_path)
#%%folder
#
foldernama='./OUTFinalProfR/'
fname=next(os.walk(foldernama))[1]





#%%pseudo data
#CDC
# hahat=glob.glob('/home/el/Data/Compil_data-Kriterres/061218-Cul-du-chien/Fit-avec-baseOUTdtrou30_rtrou4_tr5.0/twt*.csv')
# ouca='Bilb'

# fontouney=20
# listounet=[]
# for filit in hahat:

#     fin=filit.find('.csv',-4)   
#     guili=filit.find('/',-30)
#     debut=filit.find(ouca,guili)+len(ouca)
    
#     Nama=filit[debut:fin]
#     #Nama=filit[-11:-4] #Poligny
#     #Nama=filit[-12:-4]#Bilbo
#     #Nama=filit[-12:-4]#Auffargis
#     #temp=np.genfromtxt(filit, delimiter=',')
#     temp=np.genfromtxt(filit, delimiter=',',skip_header=1)
#     TWT_XP=temp[:,1]
#     Time_TWT_XP=temp[:,0]



#     temp=np.genfromtxt(filit[0:guili]+'/'+'volumes_'+ouca+Nama+'.csv',delimiter=',',skip_header=1)            
#     VOL_XP_temp=temp[:,1]
#     Time_VOL_XP=temp[:,0]
#     VOL_XP=0.001*np.interp(Time_TWT_XP,Time_VOL_XP,VOL_XP_temp)
#     listounet.append([(TWT_XP),(VOL_XP),[255, 165, 0]])

# df_params_XP=pd.DataFrame(listounet,columns=['TWT_XP','VOL_XP','Color']) 
TWT_XP=F_extractTWT('./OUTFinalProfR/OUTProfdtrou30.0_rtrou4.0_tr5.0/ts0.36_ti0.09_tr0.03_n6_alpha0.025_Ks0.4/')
VOL_XP=np.genfromtxt('./OUTFinalProfR/OUTProfdtrou30.0_rtrou4.0_tr5.0/ts0.36_ti0.09_tr0.03_n6_alpha0.025_Ks0.4'+'/Volumes_EL.csv',delimiter=',',skip_header=1)

#%%
#ouca='Poligny'
#ouca='Bilb'
#ouca='Auffar'
ouca='Tcherno'
fontouney=20
lst=[]
for ii in fname: 

    p=read_parameters(foldernama+ii+'/'+'ts0.36_ti0.09_tr0.03_n6_alpha0.025_Ks0.4/')
    paramMVG = ParamMVG(tr=p[2], ts=p[0], ti=p[1], Ks=p[5], n=p[3], alpha=p[4])
    try:
        
        temp=F_extractTWT(foldernama+ii+'/ts0.36_ti0.09_tr0.03_n6_alpha0.025_Ks0.4/')
        vol=np.genfromtxt(foldernama+ii+'/ts0.36_ti0.09_tr0.03_n6_alpha0.025_Ks0.4'+'/Volumes_EL.csv',delimiter=',',skip_header=1)
        rmseTwt=np.sqrt(np.mean((temp-TWT_XP)**2)/len(TWT_XP))/(max(TWT_XP)-min(TWT_XP))
        rmsevol=np.sqrt(np.mean(((0.001*(vol-VOL_XP))**2))/len(VOL_XP))/(max(VOL_XP)-min(VOL_XP))
            #rmsevol=0
        rmse=np.sqrt(rmseTwt**2+rmsevol**2)
        #vol=np.nan
        bibi = 0
        #rmse=np.sqrt(rmseTwt**2+rmsevol**2)
    except:
        bibi=1
        #rmseTwt=np.nan
        #rmsevol=np.nan
        rmse=np.nan
    
    lst.append([paramMVG.tr,paramMVG.ts,paramMVG.ti,paramMVG.n,paramMVG.alpha,paramMVG.Ks,bibi,rmse,(temp),(vol),float(ii[-9:-6]),float(ii[-19:-15])])
    
df_params=pd.DataFrame(lst,columns=['tr','ts','ti','n','alpha','Ks','Converged','RMSE','TWT','VOL','radius','depth']) 

df_params['alpha']=df_params['alpha']*100


#cmap = mpl.cm.jet((df_params_sorted_cut['RMSE'].values - df_params_sorted_cut['RMSE'].min())/\
#(df_params_sorted_cut['RMSE'].max() - df_params_sorted_cut['RMSE'].min()))  

  

plt.close('all')
Time_TWT_XP=[0,0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]
from matplotlib.cm import get_cmap
cmap = mpl.cm.autumn

(f2, ax)= plt.subplots(2,1,figsize=(10,10))    
k=0
ax[0].grid()
ax[1].grid()
for  (index, row) in df_params.iterrows():
    #ax[0].plot(Time_TWT_XP, row['TWT'], c=cmap(k/float(len(df_params))))
    #ax[1].plot(Time_TWT_XP, 0.001*row['VOL'], c=cmap(k/float(len(df_params))))
    ax[0].plot(Time_TWT_XP, row['TWT'], c='r')
    ax[1].plot(Time_TWT_XP, 0.001*row['VOL'], c='r')
    #k=k+1
    
# for  (index, row) in df_params_XP.iterrows():
#     ax[0].scatter(Time_TWT_XP, row['TWT_XP'],c='k',marker='*',s=100)
#     ax[1].scatter(Time_TWT_XP, row['VOL_XP'],c='k',marker='*',s=100)
   

ax[0].tick_params(axis='both', which='major', labelsize=fontouney)
ax[1].tick_params(axis='both', which='major', labelsize=fontouney)
ax[0].set_ylabel('TWT (ns)',fontsize=fontouney)
ax[1].set_ylabel('Vol (l)',fontsize=fontouney)
ax[1].set_xlabel('Experimental time (min)',fontsize=fontouney)
ax[0].set_xlabel('Experimental time (min)',fontsize=fontouney)


#f2.savefig('./plots/Synthe_Rosetta.png',format='png')


#%%







################Sensibility plot



plt.close('all')
#df_params=df_params[(df_params['tr']==0.03) & (df_params['Ks']<0.49) & (df_params['Ks']>0.07) & (df_params['n']<10.1) ]

(f1, ax)= plt.subplots(1,1,figsize=(25,15))
#cmap = mpl.cm.jet(vmin=0, vmax=1)
#norma = mpl.colors.Normalize(vmin=0, vmax=1)

#maxi=np.log10(1)
#mini=np.log10(0.05)
#maxi=2
#mini=0
#df_params['RMSE'] = df_params['RMSE'].apply(np.log10)
#norm=plt.Normalize(mini,maxi)
sc=ax.scatter(df_params['radius'],df_params['depth'],c=df_params.RMSE,cmap = 'jet')
ax.plot(4.0,30,marker='x',markersize=22)
ax.set_xlabel('R',fontsize=fontouney)
ax.set_ylabel('Depth',fontsize=fontouney)
ax.tick_params(axis='both', which='major', labelsize=fontouney)
ax.grid() 
cbar=plt.colorbar(sc)
cbar.ax.tick_params(labelsize=fontouney)         






#f1.savefig('./plots/'+ouca+'/RMSEVOLANDTWT_'+str(np.round(100*pc))+'pc_'+ouca+Nama+'.png',format='png')
#%%    
####################    
################### histogram
    # df_params_sorted=df_params.sort_values(by=['RMSE'],inplace=False)
    # df_params_sorted.reset_index(drop=True, inplace=True)
    # pc=0.1#10percent
    
    # df_params_sorted_cut=df_params_sorted[0:np.int(np.round(pc*len(df_params_sorted)))]
    # df_params_sorted_cut.reset_index(drop=True, inplace=True)
    
    # #df_params_sorted_cut=df_params_sorted[df_params_sorted['RMSE']<0.1]
    # #df_params_sorted_cut.reset_index(drop=True, inplace=True)

    # plt.close('all')
    # legendounet=['ts','n','alpha','Ks']
    # (f2, ax)= plt.subplots(2,3,figsize=(25,15))
    # kk=0
    
    # for ii in range(2):
    #     for jj in range(2):
    #         ax[ii,jj].hist(df_params_sorted_cut[legendounet[kk]], weights=np.zeros_like(df_params_sorted_cut[legendounet[kk]]) + 1. / df_params_sorted_cut[legendounet[kk]].size)
    #         ax[ii,jj].set_xlabel(legendounet[kk]+' Optim= '+str(df_params_sorted_cut[legendounet[kk]][0]),fontsize=fontouney)
    #         ax[ii,jj].set_ylabel('Rel Freq.',fontsize=fontouney)
    #         #ax[ii,jj].plot(df_params_sorted_cut[legendounet[kk]][0], 0.15, marker = '*', markersize = 10, mfc = 'k')
    #         ax[ii,jj].tick_params(axis='both', which='major', labelsize=fontouney)
    #         ax[ii,jj].grid() 
    #         ax[ii,jj].set_ylim([0, 0.3])
    #         kk=kk+1




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
    
    # cmap = mpl.cm.jet((df_params_sorted_cut['RMSE'].values - df_params_sorted_cut['RMSE'].min())/\
    # (df_params_sorted_cut['RMSE'].max() - df_params_sorted_cut['RMSE'].min()))    
    # for  color, (index, row) in zip(cmap,df_params_sorted_cut.iterrows()):
    #     ax[0,2].plot(Time_TWT_XP, row['TWT'], c=color)
    #     #ax[1,2].plot(Time_TWT_XP, row['VOL'], c=color)
    # ax[0,2].grid()
    # ax[1,2].grid()
    # ax[0,2].plot(Time_TWT_XP, TWT_XP, marker = '+', markersize = 10, mfc = 'k')
    # ax[0,2].plot(Time_TWT_XP,df_params_sorted_cut['TWT'][0], marker = '*', markersize = 10, mfc = 'k')
    # #ax[1,2].plot(Time_TWT_XP, 0.001*VOL_XP, marker = '+', markersize = 10, mfc = 'k') 
    # #ax[1,2].plot(Time_TWT_XP,df_params_sorted_cut['VOL'][0], marker = '*', markersize = 10, mfc = 'k')
    # ax[0,2].tick_params(axis='both', which='major', labelsize=fontouney)
    # ax[1,2].tick_params(axis='both', which='major', labelsize=fontouney)
    # ax[0,2].set_ylabel('TWT(ns)',fontsize=fontouney)
    # ax[0,2].set_xlabel('Exp. Time (min)',fontsize=fontouney)
    # #ax[1,2].set_ylabel('Vol. (l)',fontsize=fontouney)
    # #ax[1,2].set_xlabel('Exp. Time (min)',fontsize=fontouney)
    
    
    # f2.suptitle(str(np.round(100*pc))+'% des modèles '+ ouca+'-'+Nama+'\n Nombre de modeles: '+str(len(df_params_sorted_cut))+', RMSE: ['+str(np.round(df_params_sorted_cut['RMSE'][0],4))+';'+str(np.round(df_params_sorted_cut['RMSE'][len(df_params_sorted_cut)-1],4))+'] (ns)', fontsize=fontouney)
    # f2.savefig('./plots/'+ouca+'/Histo_'+'pc_'+str(round(100*pc))+'_'+ouca+'-'+Nama+'.png',format='png')
###################

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










