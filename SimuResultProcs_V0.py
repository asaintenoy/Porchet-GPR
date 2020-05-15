#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 22:49:07 2020

@author: el
"""
import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import datetime as dt
import os
import itertools
import io
pd.set_option('max_columns', 7)
#%% pathounet
data_path='/home/el/Codes/Porchet-GPR/OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'


X_path='/home/el/Codes/Porchet-GPR/'

hehe=os.getcwd()

if(hehe==X_path):
    print('HEHEHE on est bon')
else:
    os.chdir(X_path)
#%% Reading the folder names
fname=next(os.walk(data_path))[1]

#%% pour chaque sous folder, on lit le fichier Params
#s = "Param(a=1, b=2)"



lst=[]
for ii in fname: 

    s = io.open('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+'/Parameters').read()
    temp=(s[9:-4]).split('=')
    temp3=','.join(temp)
    temp2=temp3.split(',')
    try:
        ss = io.open('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT').read()
        bibi=0
    except:
        bibi=1
    
    lst.append([float(temp2[1]),float(temp2[3]),float(temp2[5]),float(temp2[7]),float(temp2[9]),float(temp2[11]),bibi])
    
df_params=pd.DataFrame(lst,columns=['tr','ti','ts','n','alpha','Ks','Converged'])                 

#%% seqborn
# plt.close('all')
# f1, ax1= plt.subplots(1,1,figsize=(25,15))
# sns.set(style="whitegrid")

g = sns.PairGrid(df_params, diag_sharey=False)
g.map_upper(sns.scatterplot)
g.map_lower(sns.scatterplot)
#g.map_lower(sns.kdeplot, colors="C0")
g.map_diag(sns.distplot,kde=False)


