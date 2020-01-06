import numpy as np
import pygimli as pg
from pygimli.meshtools import appendTriangleBoundary, merge2Meshes
from pygimli.mplviewer import drawMesh
from pygimli.viewer import showMesh
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
#from pygimli.mplviewer import drawMesh, drawModel
from pygimli.meshtools import interpolate
from pygimli.meshtools import nodeDataToCellData
import math

def CRIM(theta, paramMVG, paramGPRMAX):
    """Compute the relative permittivity value corresponding to the water content theta using CRIM relation"""
    return round(( math.sqrt(paramGPRMAX.eps_w)*theta+(1-paramMVG.porosity)*math.sqrt(paramGPRMAX.eps_s)+(paramMVG.porosity-theta) )**2,3)
        
def TOPP(theta):
    """Compute the permittivity value corresponding to the water content theta using the Topp relation"""
    return round(3.03+9.3*theta+146*(theta**2)-76.6*(theta**3),3)

def Rhoades(theta):
    """"Compute the electrical conductivity (S/m) corresponding to the water content theta,
    using Rhoades relation with A and B coefficients and sigma_w """
    A=1.21
    B=0.132
    sigma_w=5 
    return round((A*(theta**2)+B*theta)*sigma_w,1)

def maillage_GPRMAX(paramGPRMAX, paramMVG, mesh, mesh_pos, f_thetas, nT):    

    xmin=paramGPRMAX.xmin
    xmax=paramGPRMAX.xmax
    zmin=paramGPRMAX.zmin
    zmax=paramGPRMAX.zmax
    dx=paramGPRMAX.dx

    xreg = np.arange(xmin, xmax + dx, dx, 'float')
    zreg = np.arange(zmin, zmax + dx, dx, 'float')

    mesh2 = pg.Mesh(3)
    mesh2.createGrid(xreg, zreg)
    for c in mesh2.cells():
        c.setMarker(3)
        
    pg_pos2 = mesh2.positions()
    #On crée une matrice contenant la position des noeuds
    mesh2_pos2 = np.array((np.array(pg.x(pg_pos2)), np.array(pg.y(pg_pos2)), np.array(pg.z(pg_pos2)))).T 
    #Matrice vide de la taille du nombre de cellules
    mesh2_cells2 = np.zeros((mesh2.cellCount(), 4))

    #On rentre les cellules dans une matrice
    for i, cell in enumerate(mesh2.cells()): 
        mesh2_cells2[i] = cell.ids()
        
    mx = pg.x(mesh2.cellCenter())
    my = pg.y(mesh2.cellCenter())
    mesh_pos2=mesh2_pos2[:,0:2]
    
    xv, yv = np.meshgrid(xreg, zreg, sparse=False, indexing='ij')

    #maillage triangulaire que l'on a défini pour SWMS_2D
    maillage = mesh_pos 

    (x,z)=np.shape(maillage)
    
    theta=np.loadtxt(f_thetas) #Ouvrir le fichier
    min_theta=min(theta)
    
    eps=np.zeros(len(theta))

    for i in range(0,len(theta)):
        eps[i]=CRIM(theta[i], paramMVG, paramGPRMAX)
    
    eps_mat=np.zeros([x,int((len(eps)/x))])
    for i in range(0,nT+1) : #
        xi=i*x
        eps_mat[:,i]=eps[xi:(xi+x)]
        
    #grid_lin=np.zeros([len(mx),nT+1])
    grid_mat={}

    #fig, ax = plt.subplots(nT+1, figsize=(20, 100))
    for i in range(0, nT+1) : #
        grid=np.zeros([len(xv[:,0]), len(xv[0,:])])
        outdata=interpolate(mesh2,mesh,eps_mat[:,i], fill_value=eps[0])
        outdata2=nodeDataToCellData(mesh2,outdata)
        for j in range(0,len(xv[0,:])):
            k=j*len(xv[:,0])
            kk=len(xv[:,0])
            grid[:,j]=np.around(outdata[k:(k+kk)], decimals=2)
        grid_mat[i]=grid.T
        a=np.where(grid_mat[i]==0.0)
        grid_mat[i][a]=min(eps)
        b=np.where(grid_mat[i]<=eps[0])
        grid_mat[i][b]=eps[0]
        #drawModel(ax[i], mesh2 , outdata2)
        
    #Même interpolation pour les sigma si nécessaire
    sigma=np.zeros(len(theta))
    
    for i in range(0,len(theta)):
        #if theta[i]==min_theta:
        #    sigma[i]=0.0
        sigma[i]=Rhoades(theta[i])
    
    sigma_mat=np.zeros([x,int((len(eps)/x))])
    for i in range(0, nT+1):
        xi=i*x
        sigma_mat[:,i]=sigma[xi:(xi+x)]

    sigma_grid_mat={}
    for i in range(0, nT+1) : #
        grid=np.zeros([len(xv[:,0]), len(xv[0,:])])
        outdata=interpolate(mesh2,mesh,sigma_mat[:,i], fill_value=min(sigma))
        outdata2=nodeDataToCellData(mesh2,outdata)
        for j in range(0,len(xv[0,:])):
            k=j*len(xv[:,0])
            kk=len(xv[:,0])
            grid[:,j]=np.around(outdata[k:(k+kk)], decimals=1)
        sigma_grid_mat[i]=grid.T
        a=np.where(sigma_grid_mat[i]==sigma[0])
        sigma_grid_mat[i][a]=0.0
    
    return xv, yv, mx, my, mesh2, grid, grid_mat, eps_mat, sigma_grid_mat
