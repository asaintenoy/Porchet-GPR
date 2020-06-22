# coding: utf-8

# In[ ]:

import numpy as np
import pygimli as pg
from pygimli.meshtools import polytools as plc
def maillage_SWMS2D_EL(geometry):
    """Définition du maillage triangulaire pour SWMS_2D"""

    xmin=geometry.xmin
    xmax=geometry.xmax
    emin=geometry.emin
    emax=geometry.emax
    dtrou=geometry.dtrou
    etrou=geometry.etrou
    r=geometry.r
    dx=geometry.dx
    zaff=geometry.zaff
    eaff=geometry.eaff
    #waff=geometry.waff
    
    assert dtrou + zaff < emax

    #xtrou_reg = np.arange(xmin, r + dx, dx, 'float')
#    xtrou_reg = np.arange(xmin, r + zaff, dx, 'float')
    #etrou_reg = np.arange(etrou, emax +dx, dx, 'float')
    #efin_reg = np.arange(eaff, etrou+dx, dx, 'float')
    
    #A présent on crée une zone grâce à un polygone

    poly = pg.Mesh(2)  # empty 2d mesh
    #nStart = poly.createNode(0.0, 0.0, 0.0) # On crée un noeud de départ, on travaille en 2D donc le dernier terme vaut 0.0


    #nA = nStart # noeud de départ
    xreg=[xmin,xmin,xmin+r,xmin+r,xmax,xmax,xmin+r]
    zreg=[emin,etrou,etrou,emax,emax,emin,emin]
    nStart = poly.createNode(xreg[0], zreg[0], 0.0) # On crée un noeud de départ, on travaille en 2D donc le dernier terme vaut 0.0

    nA = nStart # noeud de départ
    for xx,zz in zip(xreg[1::],zreg[1::]): # On démarre de 1 et on se balade sur l'axe des x en créant un noeud à chaque fois
        nB = poly.createNode(xx, zz, 0.0)
        poly.createEdge(nA, nB) # On définit un côté entre le noeud précédemment créé et le nouveau
        nA = nB # On remplace le noeud de départ par le noeud nouvellement créé
    poly.createEdge(nA, nStart)

    
#=============================================================================
    c1 = plc.createCircle(pos=[xmin+r/2, etrou+1], radius=20, area=geometry.area*0.3)
    mesh1=pg.meshtools.createMesh(c1, quality=geometry.quality, area=geometry.area, smooth=geometry.smooth)
    #pg.show(mesh1, markers=True, showMesh=True)
    
    for ii in range(mesh1.nodeCount()):
        if(mesh1.node(ii)[1]>etrou):
            if(mesh1.node(ii)[0]>xmin+r):
                poly.createNode(mesh1.node(ii)[0],mesh1.node(ii)[1],0)

        elif(mesh1.node(ii)[1]<etrou):
            if(mesh1.node(ii)[0]>xmin):
                poly.createNode(mesh1.node(ii)[0],mesh1.node(ii)[1],0)




#=============================================================================

    

    mesh=pg.meshtools.createMesh(poly, quality=geometry.quality, area=geometry.area, smooth=geometry.smooth)
    #mesh=pg.meshtools.createMesh(poly, quality=34, area=0.5, smooth=geometry.smooth)
    poly.exportVTK('plc')
    mesh.exportVTK('mesh')
    pg_pos = mesh.positions()
    mesh_pos = np.array((np.array(pg.x(pg_pos)), np.array(pg.y(pg_pos)), np.array(pg.z(pg_pos)))).T #On crée une matrice contenant la position des noeuds
    mesh_cells = np.zeros((mesh.cellCount(), 3)) #Matrice vide de la taille du nombre de cellules
    for i, cell in enumerate(mesh.cells()): #On rentre les cellules das une matrice
        mesh_cells[i] = cell.ids()
    
    return mesh, pg_pos, mesh_pos, mesh_cells
