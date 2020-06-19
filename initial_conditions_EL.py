import numpy as np
from outils import comparaison_array_number,comparaison_number_number

def initial_conditions_EL(mesh_pos, geometry, paramMVG):
    """Calcule les conditions initiales du modèle hydro
    p est un tableau de coordonnees x,y, type, charge
    où le type decrit le type de noeud, charge est la charge imposée initiale
    """

    h0 = paramMVG.h0()
    
    p = np.zeros((len(mesh_pos), 4), dtype = 'float') #matrice vide
    p[:,:2] = np.array(mesh_pos[:,:2]) #coordonnées xy des noeuds du maillage

    for i in range(0,len(p)) :
        # charge au fond du trou
        if comparaison_number_number(p[i,1],geometry.etrou,geometry.tol) and p[i,0]<=(np.float64(geometry.r)-geometry.tol):
            p[i,2]=1 # type de neoud qui ne peut pas évoluer
            p[i,3]=geometry.h_eau
        # charge au fond du modèle    
        elif comparaison_number_number(p[i,1],0,geometry.tol):
            p[i,2]=0 #type free drainage?
            p[i,3]=h0
        # partout ailleurs    
        else :
            p[i,3]=h0
                
    #Interpolation lineaire pour la valeur de charge au bord du trou sous l'eau
    #On considère que la charge est nulle à la surface de l'eau
    rows = np.array(np.where((p[:,1] <= np.float64(geometry.etrou) + np.float64(geometry.h_eau) - geometry.tol) & (p[:,1] >= np.float64(geometry.etrou) + geometry.tol )\
                             & (comparaison_array_number(p[:,0],geometry.r,geometry.tol)))).T
    
    for i in range(1,len(rows)):
        p[rows[i,0],2] = 1
        p[rows[i,0],3] = np.float64(geometry.h_eau)-(p[rows[i,0],1]-np.float64(geometry.etrou))
     
    return p

