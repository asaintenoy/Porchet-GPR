import h5py
import math
import numpy as np
from scipy.stats import linregress

def picking(filename, A_tab, nT, geometry) :
    """ Search for TWT of the wave going around the bulb in each traces.
    The minimum TWT is computed using the geometry and min relative permittivity.
    On each trace, it searches the max amplitude TWT. Then the min TWT just before and the min TWT just after.
    It fits the three curves TWT(itrace) with a power law and take the one with the best fit. 
    """
    #TODO:ajuster les tailles de fenetre... C'est un peu bizarre pour le moment.
    
    f = h5py.File(filename, 'r')
    path = '/rxs/rx1/'
    samples = f.attrs['Iterations']
    dt = f.attrs['dt']*1e9 # en ns
    dx = 1
    data = np.ones((samples, nT+1))
    t_max = np.zeros(nT+1)
    t_min1 = np.zeros(nT+1)
    t_min2 = np.zeros(nT+1)
    tt=np.zeros(nT+1)
    tt_min1 = np.zeros(nT+1)
    tt_min2 = np.zeros(nT+1)
    tps_max = np.zeros(nT+1)
    tps_min1 = np.zeros(nT+1)
    tps_min2 = np.zeros(nT+1)

    # Estimate minimun TWT for the first reflexion
    eps_init = A_tab[0][:,2].min()
    #h = math.sqrt(((geometry.dtrou-geometry.h_eau)*0.01)**2 + ((geometry.d_emet+geometry.d_recept)/2)**2) #Pourquoi la distance trou-antenne n'est pas mise en m???
    h = math.sqrt(((geometry.dtrou-geometry.h_eau)*0.01)**2 + ((paramGPRMAX.d_emet+paramGPRMAX.d_recept)*0.01/2)**2) # en m
    v_init=0.3/(math.sqrt(eps_init)) # en m/ns
    t_init = (2*h)/v_init # en ns
    print(t_init)
    itmin0 = t_init/dt
    fenetre = 2 #taille de la fenetre de recherche de reflexion max ou min (en ns)
    ifenetre = int(fenetre/dt)

    #Picking onde du bas (onde qui fait le tour du bulbe)

    for itrace in range(0,nT+1):
        data[:,itrace] = f['%s%s' % (path, 'Ez')][:,itrace]
        if itrace==0 :
            itmin = itmin0
        else :
            itmin = tt[itrace-1]
            
        t_max[itrace] = np.max(data[int(itmin):int(itmin+ifenetre*3.5/(itrace+1)),itrace])
        tt[itrace] = np.where(data[:,itrace]==t_max[itrace])[0]
        tps_max[itrace] = tt[itrace]*dt
        t_min1[itrace]= np.min(data[int(tt[itrace]-ifenetre*(itrace+1)/4):int(tt[itrace]),itrace])
        tt_min1[itrace] = np.where(data[:,itrace]==t_min1[itrace])[0]
        tps_min1[itrace] = tt_min1[itrace]*dt
        t_min2[itrace] = np.min(data[int(tt[itrace]):int(tt[itrace]+ifenetre),itrace])
        tt_min2[itrace] = np.where(data[:,itrace]==t_min2[itrace])[0]
        tps_min2[itrace] = tt_min2[itrace]*dt

    tps_max = tps_max - tps_max[0]
    tps_min1 = tps_min1 - tps_min1[0]
    tps_min2 = tps_min2 - tps_min2[0]
    
    # Ajusting a power law
    def power_law(t,a,b) :
        return np.exp(a*np.log(t)+b)

    def rmse(predictions, targets):
        return np.sqrt(((predictions-targets) ** 2).mean())

    #On calcule la regression linéaire sur log(tps_bas) et log(tps_bas2)
    Sol_max=linregress(np.log(temps),np.log(tps_max[1:]))
    Sol_min1=linregress(np.log(temps),np.log(tps_min1[1:]))
    Sol_min2=linregress(np.log(temps),np.log(tps_min2[1:]))

    #On calcule les temps d'arrivée du modèle correspondant à la regression
    Model_TWT_max=power_law(temps,Sol_max[0],Sol_max[1])
    Model_TWT_min1=power_law(temps,Sol_min1[0],Sol_min1[1])
    Model_TWT_min2=power_law(temps,Sol_min2[0],Sol_min2[1])

    #On calcule la RMSE entre TWT de la fonction puissance et TWT trace  
    rmse_max=rmse(Model_TWT_max,tps_max[1:])
    rmse_min1=rmse(Model_TWT_min1,tps_min1[1:])
    rmse_min2=rmse(Model_TWT_min2,tps_min2[1:])

    #On compare les deux rmse et l'on selectionne l'un ou l'autre des deux pickings
    if rmse_max>rmse_min1:
        if rmse_min1 > rmse_min2:
            twt_fin = tps_min2
            cas = 'min2'
        else:
            twt_fin = tps_min1
            cas = 'min1'
    else:
        if rmse_max > rmse_min2:
            twt_fin = tps_min2
            cas = 'min2'
        else:
            twt_fin = tps_max
            cas = 'max'

    return cas, twt_fin
