#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:42:18 2020

@author: el
"""



import os
import matplotlib.pyplot as plt
from modelisation import run
from param_acquisition import Geometry, ParamMVG, ParamGPRMAX
from picking_radargramme import picking

def Forward(geometry,paramMVG,paramGPRMAX,temps,tmax_SWMS2D):
#%%
    folderout=run(geometry=geometry,paramMVG=paramMVG,paramGPRMAX=paramGPRMAX,temps=temps,tmax_SWMS2D=tmax_SWMS2D)
    try:
        cas,dt,itmin0,ifenetre,tps_min1,tps_min1_0,tps_min2,tps_min2_0,tps_max,tps_max0,TWT\
            =picking(os.path.join(folderout, 'radargram__merged.out'), len(temps), geometry, paramMVG, paramGPRMAX, temps)
    except:
        print('Probably already simulated')
    return TWT
