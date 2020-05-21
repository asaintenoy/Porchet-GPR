#import math
#from scipy.interpolate import griddata
import numpy as np
import os
import time
retval = os.getcwd()
#import matplotlib.pyplot as plt
#import sys, h5py, binascii
#import tempfile
#import pygimli as pg
#from pygimli.meshtools import appendTriangleBoundary, merge2Meshes, mesh
import warnings

#Pour aider pyflakes à analyser il vaut mieux importer les fonctions explicitement
from maillage_SWMS2D_EL import maillage_SWMS2D_EL
from initial_conditions import initial_conditions
from ecriture_fichiers_SWMS2D import ecriture_Selector_in, ecriture_Grid_in 
from ecriture_fichiers_GPRMAX import ecriture_fichiers_GPRMAX
from maillage_GPRMAX import maillage_GPRMAX
#from outils import showQuality
#from pygimli.meshtools import quality
#from picking_radargramme import picking

from param_acquisition import longueur_d_onde 
from joblib import Parallel, delayed


def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

#Def des paramètres de géométrie du modéle
# class Geometry :
#     def __repr__(self):
#         return "dtrou{dtrou}_rtrou{r}_tr{h_eau}".format(dtrou=self.dtrou,r=self.r,h_eau=self.h_eau)        

# # Definition des paramètres MVG
# class ParamMVG(namedtuple("ParamMVG", ["ts", "ti", "tr", "n", "alpha", "Ks"])):
               
#     #Recalcul des h0 à partir des param MVG
#     def h0(self):
#         return -(1/self.alpha)*(((self.ti-self.tr)/(self.ts-self.tr))**(self.n/(1-self.n))-1)**(1/self.n)
#     def __repr__(self):
#         return "ts{ts}_ti{ti}_tr{tr}_n{n}_alpha{alpha}_Ks{Ks}".format(ts=self.ts,ti=self.ti,tr=self.tr,n=self.n,alpha=self.alpha,Ks=self.Ks)        

# # Definition des param gprMax
# class ParamGPRMAX :
#     pass

# def longueur_d_onde(theta, freq, paramMVG, paramGPRMAX):
#     """Compute the wavelength for a given water content and frequency (Hz)"""
#     eps = CRIM(theta,paramMVG,paramGPRMAX)
#     vitesse = 0.3 / math.sqrt(eps) # m/ns
#     return vitesse*(10**9)/freq

def run(geometry,paramMVG,paramGPRMAX,temps,tmax_SWMS2D):
    # try: c'est pour pouvoir prendre la main si ça foire.
    #try:

    myDirName = "OUTTEST"+repr(geometry)+"/"+repr(paramMVG)
    nom='radargram'
    filename = nom + '__merged.out'

    dir = os.getcwd() #répertoire où on a lancer le script

    if os.path.isfile(myDirName+"/"+filename) :
        print('already simulated')
        return
    
    os.makedirs(myDirName,exist_ok=True)
    os.chdir(myDirName)

    fparam=open("Parameters","w")
    fparam.write("""ParamMVG(tr={tr},ti={ti},ts={ts},n={n},alpha={alpha},Ks={Ks})\n""".format(ts=paramMVG.ts,ti=paramMVG.ti,tr=paramMVG.tr,n=paramMVG.n,alpha=paramMVG.alpha,Ks=paramMVG.Ks))
    fparam.close()
                 
    os.system("cp "+dir+"/gprMax .")
    os.system("cp "+dir+"/gprMaxMerge .")
    # ce serait mieux de copier l'executable H2D dans le dossier où l'on travaille mais il est gros!
    #os.system("cp "+dir+"/HD2/H2D .")
    
    os.makedirs("SWMS_2D.IN",exist_ok=True)
    os.makedirs("SWMS_2D.OUT",exist_ok=True)
    
    #Définition du maillage triangulaire pour SWMS2D
    #[mesh, pg_pos, mesh_pos, mesh_cells]=maillage_SWMS2D(geometry)
    [mesh, pg_pos, mesh_pos, mesh_cells]=maillage_SWMS2D_EL(geometry)
    #from pygimli.meshtools import mesh
    #from pygimli.mplviewer import drawMesh
    
    #from pygimli.viewer import showMesh
    #showMesh(mesh)
    
    #matplotlib.use('Agg')
    #figi=showQuality(mesh, quality(mesh))
    #figi.savefig('mesh.png',format='png')

    #Calcul des charges initiales en chaque noeud du maillage
    p=initial_conditions(mesh_pos, geometry, paramMVG)

    #Création des fichiers .in nécessaires pour SMWS_2D
    ecriture_Selector_in(mesh, paramMVG, temps, p)
    ecriture_Grid_in(mesh, p)
    os.system("mv Grid.in SWMS_2D.IN")
    os.system("mv Selector.in SWMS_2D.IN")
            
    #Lancement SWMS_2D
    #os.system("/home/clemence/Porchet-GPR/source/HD2/H2D")
    start_running = time.time()
    end_running=[]
    error=os.system("timeout {} {}/HD2/H2D".format(tmax_SWMS2D,dir))
    if error : #reagira seulement si error est différent de 0
        print("HD2 fut trop long")
        os.chdir(dir)
        
    end_running.append(time.time()-start_running)

    #os.popen("rm -rf *.in")
    print('Running Hydrus time:'+str(end_running))   
     

    #Fichier contenant les thetas
    f_thetas = "SWMS_2D.OUT/th.out" 

    #Number of traces to be calculated
    nT=len(temps)

    #Interpolation du maillage triangulaire sur une grille rectangulaire pour gprMax
    #On crée un maillage rectangulaire avec les dimensions du modèle
    [xv, yv, mx, my, mesh2, grid, grid_mat, eps_mat, sigma_grid_mat] = maillage_GPRMAX(paramGPRMAX, paramMVG, mesh, mesh_pos[:,:2], f_thetas, nT)
   
    # Lancement gprMax

    # pas de calcul spatial le plus grossier possible (m)
    dl = longueur_d_onde(paramMVG.ts,paramGPRMAX.freq_max,paramMVG,paramGPRMAX)/paramGPRMAX.spatial_step

    materiaux = {}
    A_tab={}

    def materiau(x, sigma): # TODO: A regarder car je n'y comprend rien!
        if x in materiaux: # si un materiau existe déjà pour x, il le renvoie, donc il n'en crée pas un nouveau.
            return materiaux[x] 
        valeur1 = "sand{}".format(len(materiaux))
        valeur2 = sigma
        materiaux[x] = valeur1, valeur2
        return valeur1, valeur2

    xreg = np.arange(paramGPRMAX.xmin, paramGPRMAX.xmax + paramGPRMAX.dx, paramGPRMAX.dx, 'float')
    zreg = np.arange(paramGPRMAX.zmin, paramGPRMAX.zmax + paramGPRMAX.dx, paramGPRMAX.dx, 'float')
    end_running=[]
    start_running = time.time()
    def rungprmax(ite,xreg,zreg,grid_mat,sigma_grid_mat,xv,yv,nom, paramMVG, paramGPRMAX, geometry, dl,materiaux):
        from ecriture_fichiers_GPRMAX import ecriture_fichiers_GPRMAX

        for j in range(0,len(zreg)):
            for k in range(0,len(xreg)):
                materiau(grid_mat[ite][j,k], sigma_grid_mat[ite][j,k])
        
        A = ecriture_fichiers_GPRMAX(xv.T*0.01, yv.T*0.01, grid_mat[ite], ite, nom, paramMVG, paramGPRMAX, geometry, dl, materiaux) 
        fichier=nom+'_'+str(ite+1)+'.in'
        
        command="./gprMax "+fichier
        os.popen(command).readlines()
        
    Parallel(n_jobs=8)(delayed(rungprmax)(ite,xreg,zreg,grid_mat,sigma_grid_mat,xv,yv,nom, paramMVG, paramGPRMAX, geometry, dl,materiaux) for ite in range(0,nT+1))





# =============================================================================
#     for i in range(0,nT+1):
#         #start_ecriture = time.time()
#         for j in range(0,len(zreg)):
#             for k in range(0,len(xreg)):
#                 materiau(grid_mat[i][j,k], sigma_grid_mat[i][j,k])
# 
#                 
#         # Writing a file.in for running GPRMAX
#         A = ecriture_fichiers_GPRMAX(xv.T*0.01, yv.T*0.01, grid_mat[i], i, nom, paramMVG, paramGPRMAX, geometry, dl, materiaux) 
#        
#         # plt.close('all')
#         # fig = plt.figure()
#         # ax = fig.add_subplot(111)
#         # gni=ax.scatter(A[:,0],A[:,1],s=30,c=A[:,2])
#         # cbar=fig.colorbar(gni,label='Eps' )
#         # cbar.minorticks_on()
#         # fig.savefig('Gnard'+str(i)+'.png',format='png')
#         
#         #A_tab[i]=A
#         #end_ecriture.append(time.time()-start_ecriture)
#         #Lancement calcul gprMax
#         #fichier=nom+'_'+str(i+1)+'.in'
#         fichier=nom+'_'+str(i+1)+'.in'
#         
#         command="./gprMax "+fichier
#         os.popen(command).readlines()
# =============================================================================
        
        

    # Concatenate all nT traces    
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax2 = ax1.twinx()
    # ax1.plot(list(range(0,nT+1)),end_running,c='r')
    # ax2.plot(list(range(0,nT+1)),end_ecriture,c='b')
    # ax1.grid()
    # ax2.grid()
    command2="./gprMaxMerge "+ nom + "_" 
    os.popen(command2).readlines()
    end_running.append(time.time()-start_running)

    #os.popen("rm -rf *.in")
    print('Running GPR time:'+str(end_running))
    for i in range(0,nT+1) :
        os.popen("rm -rf "+ nom + "_" + str(i+1) + ".out")

    # Picking of the nT TWTs
    #cas, TWT = picking(filename, A_tab, nT, geometry)

    os.chdir(dir)
    
    return myDirName
    
#run(0.3,10,0.024,0.097,0.401,) #Ks, n, alpha, tr, ts, ti
