# coding: utf-8

#TODO: Indiquer le path de H2D en tête pour quand on change d'ordi...

import os

from modelisation import Geometry, ParamMVG, ParamGPRMAX, run

geometry=Geometry()

#Domaine de calcul (en cm)
# largeur
geometry.xmin=0 
geometry.xmax=40
# hauteur (elevation)
geometry.emin=0
geometry.emax = 80
# profondeur du trou en cm
geometry.dtrou = 30
# elevation du fond du trou
geometry.etrou = geometry.emax - geometry.dtrou
 # rayon du trou en cm
geometry.r=4
# hauteur d'eau imposée au fond du trou en cm
geometry.h_eau=10.0 
# pas de la maille en cm
geometry.dx = 0.1
# profondeur sous le trou (cm) jusqu'où on souhaite un maillage affiné. 
geometry.zaff= 20
#largeur horizontal de la zone affinée (cm)
geometry.waff=20
# elevation de l'affinage
geometry.eaff=geometry.etrou-geometry.zaff 
# contrainte d'angle min pour mesh 
geometry.quality=33
# maximum triangle size  (m*²)
geometry.area=5
# tupple for mesh generation 
geometry.smooth=[1,5]


#Temps d'infiltration où à lieu le calcul de chaque trace  (minutes) 
#temps=[1.00, 2.00, 3.00, 4.00, 5.00, 6.00]
temps=[0.17, 0.33, 0.50, 0.67, 0.83, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00]

# Temps max de calcul SWMS2D au delà duquel on arrète le calcul (secondes)
tmax_SWMS2D = 600

# Definition des paramètres MVG
paramMVG=ParamMVG()

#Teneur en eau résiduelle
paramMVG.tr = 0.03
#Teneur en eau à saturation
#paramMVG.ts=np.arange(0.36, 0.46, 0.02, 'float')
paramMVG.ts = 0.36
#Teneur en eau initiale
#paramMVG.ti=np.arange(0.05, 0.13, 0.02, 'float')
paramMVG.ti = 0.05
#Perméabilité à saturation
#paramMVG.Ks=np.arange(0.01, 0.56, 0.05, 'float')
paramMVG.Ks = 0.1
#param fitting retention n
#paramMVG.n=np.arange(1.5, 10.25, 0.25, 'float')
paramMVG.n = 4
#param fitting retention alpha
#paramMVG.alpha=np.arange(0.01, 0.11, 0.01, 'float')
paramMVG.alpha = 0.01
#Porosité
paramMVG.porosity = paramMVG.ts

# Definition des param gprMax
paramGPRMAX=ParamGPRMAX()
# Domaine de calcul (cm)
paramGPRMAX.xmin = geometry.xmin
paramGPRMAX.xmax = geometry.xmax
paramGPRMAX.zmin = geometry.emin
paramGPRMAX.zmax = geometry.emax
# Taille des mailles (cm)
paramGPRMAX.dx = 1.0 
# Electrical conductivity of the medium
paramGPRMAX.sigma=0.0000
# Relative dielectric permittivity of water
paramGPRMAX.eps_w=80.1
# Relative dielectric permittivity of PVC
paramGPRMAX.eps_pvc=3
# Relative dielectric permittivity of pure silice
paramGPRMAX.eps_s=2.5
# Ricker signal central frequency (Hz)
paramGPRMAX.wave_freq = 1000e6
# Frequence max du signal EM (Hz)
paramGPRMAX.freq_max = 2.8 * paramGPRMAX.wave_freq
# Distance between hole middle and source (m)
paramGPRMAX.d_emet = 0.18
# Distance between hole middle and receiving antenna (m)
paramGPRMAX.d_recept = 0.22
# param qui raffine le pas spatial (par défaut 10 d'après doc gprmax)
paramGPRMAX.spatial_step = 5
# Trace time window (ns)
paramGPRMAX.time = 30e-9
#time_step_stability_factor (pas utilisé pour le moment...)
paramGPRMAX.fac_dt = 0.2 

myDirName = "OUT"+repr(geometry)+"/"+repr(paramMVG)

#myDirName = tempfile.mkdtemp() #Création d'un dossier temporaire pour stocker les données et fichiers
#print(myDirName)


run(geometry=geometry,paramMVG=paramMVG,paramGPRMAX=paramGPRMAX,temps=temps,tmax_SWMS2D=tmax_SWMS2D,myDirName=myDirName)
