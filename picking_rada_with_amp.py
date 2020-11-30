import h5py
import math
import numpy as np
from scipy.stats import linregress
from maillage_GPRMAX import CRIM

def picked_amp(filename, nT, geometry, paramMVG, paramGPRMAX, temps):
    """ Get the data amplitude amp in filename at given times TWT read in TWT_EL.csv. 
    """
    #display(filename)
    f = h5py.File(filename, 'r')
    path = '/rxs/rx1/'
    samples = f.attrs['Iterations']
    dt = f.attrs['dt']*1e9 # en ns
    data = np.ones((samples, nT+1))

    tt=np.zeros(nT+1)
    tt_min1 = np.zeros(nT+1)
    tt_min2 = np.zeros(nT+1)
    tps_max = np.zeros(nT+1)
    tps_min1 = np.zeros(nT+1)
    tps_min2 = np.zeros(nT+1)
    amp_max = np.zeros(nT+1)
    amp_min1 = np.zeros(nT+1)
    amp_min2 = np.zeros(nT+1)
    
   
   # Estimate minimun TWT for the first reflexion
    eps_init = CRIM(paramMVG.ti, paramMVG, paramGPRMAX)
    #print('eps_init',eps_init)
    #h = math.sqrt(((geometry.dtrou-geometry.h_eau)*0.01)**2 +\
    h = math.sqrt((geometry.dtrou*0.01)**2 +\
                  ((paramGPRMAX.d_emet+paramGPRMAX.d_recept)/2)**2) # en m
    #print('h',h)
    v_init=0.3/(math.sqrt(eps_init)) # en m/ns
    t_init = (2*h)/v_init # en ns
    it_init = int(t_init/dt)
    
    #recherche du delta_init sur la trace 0 (temps d'arrivée du premier max sur l'onde directe)
    data[:,0] = f['%s%s' % (path, 'Ez')][:,0]
    delta_init = np.argmax(data[:it_init,0])
    
    itmin0 = it_init + delta_init
    
    periode = (1/paramGPRMAX.wave_freq)*1e9 #ns 
    #print('periode',periode)
    fenetre = 2*periode #taille de la fenetre de recherche de reflexion max ou min (en ns)
    ifenetre = int(fenetre/dt)

    for itrace in range(0,nT+1):
        data[:,itrace] = f['%s%s' % (path, 'Ez')][:,itrace]
        #print(data)
        if itrace==0 :
            itmin = itmin0
        else :
            itmin = tt[itrace-1]
            
        tt[itrace] = itmin + np.argmax(data[int(itmin):int(itmin+ifenetre),itrace])
        tps_max[itrace] = tt[itrace]*dt
        amp_max[itrace] = data[np.int(tt[itrace]),itrace]
                
        tt_min1[itrace] = tt[itrace] - ifenetre + np.argmin(data[int(tt[itrace]-ifenetre):int(tt[itrace]),itrace])
        tps_min1[itrace] = tt_min1[itrace]*dt
        amp_min1[itrace] = data[np.int(tt_min1[itrace]),itrace]
        
        tt_min2[itrace] = tt[itrace] + np.argmin(data[int(tt[itrace]):int(tt[itrace]+ifenetre),itrace])
        tps_min2[itrace] = tt_min2[itrace]*dt
        amp_min2[itrace] = data[np.int(tt_min2[itrace]),itrace]
        
#    tps_max0 = tps_max[0]
#    tps_min1_0 = tps_min1[0]
#    tps_min2_0 = tps_min2[0]
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
            amp_fin = amp_min2
 #           cas = 'min2'
        else:
            twt_fin = tps_min1
            amp_fin = amp_min1
#            cas = 'min1'
    else:
        if rmse_max > rmse_min2:
            twt_fin = tps_min2
            amp_fin = amp_min2
#            cas = 'min2'
        else:
            twt_fin = tps_max
            amp_fin = amp_max
#            cas = 'max'

    return twt_fin, amp_fin
