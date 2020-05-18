#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 09:10:58 2020

@author: el
"""
import itertools
import numpy as np
import h5py
import math
import os
import matplotlib.pyplot as plt
from scipy.stats import linregress
import io
import pandas as pd
import sys
#%%

data_path='/home/el/Codes/Porchet-GPR/OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'


X_path='/home/el/Codes/Porchet-GPR/'

hehe=os.getcwd()

if(hehe==X_path):
    print('HEHEHE on est bon')
else:
    os.chdir(X_path)
#%% Reading the folder names
fname=next(os.walk(data_path))[1]
#%%

lst=[]
blib=pd.DataFrame(columns=['TWT(ns)','Volumes'])
for ii in fname:
        
    filename = 'TWT'
    if (os.path.isfile('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+"/"+filename)):
        mots=[]
        twts=np.zeros(11)
        fTWT=open('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+"/"+filename,"r")
        i=0
        for ligne in fTWT:
            ligne = ligne.rstrip(']\n')
            mots = ligne.split(" ")
            for mot in mots:
                if (mot != '' and mot != '[0.' and mot != '[' and mot != '0.'):
                    twts[i] = float(mot)
                    i=i+1
        fTWT.close()
        twts=insert(twts,0,0)#jte jure ...
        blib['TWT(ns)']=twts
        blib['TWT(ns)'].to_csv('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+'/TWT_EL.csv', sep=',', encoding='utf-8',index=None,header=True)
     
    filename = 'Volumes'
    if (os.path.isfile('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+"/"+filename)):   
        volumes=np.zeros(11)
        fVOL=open('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+"/Volumes","r")
        i=0
        for ligne in fVOL:
              ligne = ligne.rstrip(']\n')
              mots = ligne.split(" ")
              for mot in mots:
                  if (mot != '' and mot != '[0.' and mot != '[' and mot != '0.'):
                      volumes[i] = float(mot)
                      i=i+1
        fVOL.close()  
        volumes=insert(volumes,0,0)
        blib['Volumes']=volumes
        blib['Volumes'].to_csv('./OUTdtrou30_rtrou4_tr5.0/OUTdtrou30_rtrou4_tr5.0/'+ii+'/Volumes_EL.csv', sep=',', encoding='utf-8',index=None,header=True)
        
        

        