import h5py
import math
import numpy as np
import matplotlib.pyplot as plt


# Ouverture du fichier avec le radargramme

filename = '/home/clemence/Albane-Porchet-GPR-mod/Test5/Forward279__merged.out'
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
