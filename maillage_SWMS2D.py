# coding: utf-8

# In[ ]:

import numpy as np
import pygimli as pg

def maillage_SWMS2D(geometry):
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

    xtrou_reg = np.arange(xmin, r + dx, dx, 'float')
#    xtrou_reg = np.arange(xmin, r + zaff, dx, 'float')
    etrou_reg = np.arange(etrou, emax +dx, dx, 'float')
    efin_reg = np.arange(eaff, etrou+dx, dx, 'float')
    
    #A présent on crée une zone grâce à un polygone

    poly = pg.Mesh(2)  # empty 2d mesh
    nStart = poly.createNode(0.0, 0.0, 0.0) # On crée un noeud de départ, on travaille en 2D donc le dernier terme vaut 0.0

    nA = nStart # noeud de départ

    for e in efin_reg[0:len(efin_reg)]: # On démarre de 1 et on se balade sur l'axe des x en créant un noeud à chaque fois
        nB = poly.createNode(0, e, 0.0)
        poly.createEdge(nA, nB) # On définit un côté entre le noeud précédemment créé et le nouveau
        nA = nB # On remplace le noeud de départ par le noeud nouvellement créé
    
    #Définir les noeuds au fond du trou    
    for x in xtrou_reg[1:len(xtrou_reg)]:
        nB = poly.createNode(x, etrou, 0.0)
        poly.createEdge(nA, nB) #On définit un côté entre le noeud précédemment crée et le nouveau
        nA = nB #On remplace le noeud de départ par le noeud nouvellement crée

    #Définir les noeuds le long du trou
    for f in etrou_reg[1:len(etrou_reg)]: # On démarre du haut et on se balade sur l'axe des z en créant un noeud à chaque fois
        nB = poly.createNode(r, f, 0.0)
        poly.createEdge(nA, nB) # On définit un côté entre le noeud précédemment créé et le nouveau
        nA = nB # On remplace le noeud de départ par le noeud nouvellement crée

    nC=nB
    nD = poly.createNode(xmax, emax, 0.0)
    poly.createEdge(nC, nD)
    nE = poly.createNode(xmax, emin, 0.0)
    poly.createEdge(nD, nE)
    poly.createEdge(nE, nStart) #On ferme le polygone!

    mesh=pg.meshtools.createMesh(poly, quality=geometry.quality, area=geometry.area, smooth=geometry.smooth)

    pg_pos = mesh.positions()
    mesh_pos = np.array((np.array(pg.x(pg_pos)), np.array(pg.y(pg_pos)), np.array(pg.z(pg_pos)))).T #On crée une matrice contenant la position des noeuds
    mesh_cells = np.zeros((mesh.cellCount(), 3)) #Matrice vide de la taille du nombre de cellules
    for i, cell in enumerate(mesh.cells()): #On rentre les cellules das une matrice
        mesh_cells[i] = cell.ids()
    
    return mesh, pg_pos, mesh_pos, mesh_cells
