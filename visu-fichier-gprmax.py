import h5py
import math
import numpy as np
import matplotlib.pyplot as plt


# Ouverture du fichier avec le radargramme

filename = '/home/el/Codes/Porchet-GPR/OUTdtrou30_rtrou4_tr5.0/ts0.4_ti0.1_tr0.03_n5_alpha0.03_Ks0.2/radargram__merged.out'
f = h5py.File(filename, 'r')
# Lecture des données
path = '/rxs/rx1/'
data = f['%s%s' % (path, 'Ez')][:,:]

# Affichage du radargramme

plt.figure(figsize=(20,15))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.imshow(data[int(itmin0):int(itmin0+7*ifenetre),:],aspect=0.01)
plt.imshow(data[:,:],aspect=0.01)
