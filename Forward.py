#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:42:18 2020

@author: el
"""



import os
from modelisation import run
#from param_acquisition import Geometry, ParamMVG, ParamGPRMAX
from picking_radargramme import picking
from F_extract_volumes import F_extract_volumes
import pandas as pd
import numpy as np

def Forward(geometry,paramMVG,paramGPRMAX,temps,tmax_SWMS2D):
#%%
    Res=pd.DataFrame(columns=['TWT(ns)','Volumes(ml)'])#pas le plus rapide sans doute...

    folderout=run(geometry=geometry,paramMVG=paramMVG,paramGPRMAX=paramGPRMAX,temps=temps,tmax_SWMS2D=tmax_SWMS2D)
    try:
        #cas,dt,itmin0,ifenetre,tps_min1,tps_min1_0,tps_min2,tps_min2_0,tps_max,tps_max0,TWT\
        TWT=picking(os.path.join(folderout, 'radargram__merged.out'), len(temps), geometry, paramMVG, paramGPRMAX, temps)
    except:
        print('Probably already simulated')
        TWT=np.nan
        
    Vol=F_extract_volumes(folderout,temps)
      
#pas forcement bien placé, a voir.
    Res['TWT(ns)']=TWT
    Res['Volumes(ml)']=Vol
    Res['TWT(ns)'].to_csv(folderout + '/TWT_EL.csv',sep=',', encoding='utf-8', index=None, header=True)
    Res['Volumes(ml)'].to_csv(folderout + '/Volumes_EL.csv', sep=',', encoding='utf-8', index=None, header=True)
    # os.popen('rm -rf ' + folderout + '/*.in')
    # os.popen('rm -rf ' + folderout + '/SWMS_2D.IN/*.in')
    # os.popen('rm -rf ' + folderout + '/gprMax')
    # os.popen('rm -rf ' +folderout + '/gprMaxMerge')
    # os.popen('bzip2 -9 '+folderout+ '/radargram__merged.out SWMS_2D.OUT/th.out')
    # os.popen('rm -rf ' + folderout + '/SWMS_2D.OUT/*.out')
      
    return TWT,Vol 
