#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 14:24:20 2020

@author: el
"""
import matplotlib.pyplot as plt
import numpy as np
import io
import h5py
import os

import pygimli as pg
from pygimli.meshtools import polytools as plc
from pygimli.meshtools import quality
from param_acquisition import Geometry, ParamMVG, ParamGPRMAX

def showQuality(mesh, qualities):
    fig, axes = plt.subplots(1, 2)
    axes[1].hist(qualities, color="grey")
    pg.show(mesh, qualities, ax=axes[0], cMin=0.5, cMax=1, hold=True,
            logScale=False, label="Mesh quality", cmap="RdYlGn", showMesh=True)
    s = "Min: %.2f, Mean: %.2f, Max: %.2f" % (
        np.min(qualities), np.mean(qualities), np.max(qualities))
    axes[1].set_title(s)
    axes[1].set_xlabel("Mesh quality")
    axes[1].set_ylabel("Frequency")
    axes[1].set_xlim(0, 1)

    # Figure resizing according to mesh dimesions
    x = mesh.xmax() - mesh.xmin()
    y = mesh.ymax() - mesh.ymin()
    width, height = fig.get_size_inches()
    fig.set_figheight(height * 1.3 * (y / x))
    fig.tight_layout()
    return fig

####### Lire les fichiers Parameters ######
def read_parameters(filepath):
    filename = 'Parameters'
    s = io.open(os.path.join(filepath,filename)).read()    
    p = eval(s)
    return p
 
######## Lire les fichiers TWT #############
def read_TWT(filepath):           
    filename = 'TWT'
    mots=[]
    twts=np.zeros(nT)
    fTWT=open(filepath + filename,"r")
    i=0
    for ligne in fTWT:
        ligne = ligne.rstrip(']\n')
        mots = ligne.split(" ")
        for mot in mots:
            if (mot != '' and mot != '[0.' and mot != '[' and mot != '0.'):
                twts[i] = float(mot)
                i=i+1
    return twts
    fTWT.close()
                    
### Lire les fichiers Volumes ##############
def read_volumes(filepath):
    volumes=np.zeros(nT)
    filename = 'Volumes'
    fVOL=open(filepath + filename,"r")
    i=0
    for ligne in fVOL:
        ligne = ligne.rstrip(']\n')
        mots = ligne.split(" ")
        for mot in mots:
            if (mot != '' and mot != '[0.' and mot != '[' and mot != '0.'):
                volumes[i] = float(mot)
                i=i+1
    return volumes
    fVOL.close()

### Plotter un radargramme
def rada_plot(filepath):
    filename = 'radargram__merged.out'
    f = h5py.File(filepath + filename, 'r')
    path = '/rxs/rx1/'
    data = f['%s%s' % (path, 'Ez')][:,:]
    plt.figure(figsize=(10,15))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.imshow(data[:,:],aspect=0.005)
