import numpy as np
import math
from outils import comparaison_array_number,comparaison_number_number

def ecriture_fichiers_GPRMAX_EL(X, Y, grid_z0, trace_number, nom, paramMVG, paramGPRMAX, geometry, dl, materiaux) :

    #attention, veut les param en m!!!!
    etrou = np.float64(geometry.etrou*0.01)
    h_eau = np.float64(geometry.h_eau*0.01)
    radius = np.float64(geometry.r*0.01)
    Ks = paramMVG.Ks
    sigma = paramGPRMAX.sigma
    nmedia=len(materiaux)+2
    dx=np.float64(abs((X)[0,1]-(X)[0,0]))
    dy=np.float64(abs((Y)[1,0]-(Y)[0,0]))

    if(dl>dx):
        dl=dx
        
    fgrid=open(nom+'_'+str(trace_number+1)+'.in',"w")
    fgrid.write("""------------------------------------------------\n""")
    #fgrid.write("""#number_of_media: {}\n""".format(nmedia+4))
    fgrid.write("""#domain: {:.5f} {:.5f} {:.5f}\n""".format(2*max(X[0,:])+dx, 2*max(Y[:,0]), dl))
    fgrid.write("""#dx_dy_dz: {:.5f} {:.5f} {:.5f}\n""".format(dl, dl, dl))
    fgrid.write("""#time_window: {}\n""".format(paramGPRMAX.time))
    #fgrid.write("""#time_step_stability_factor: {}\n""".format(fac_dt))
    #fgrid.write("""#abc_type: pml\n""")
    #fgrid.write("""#pml_cells: {}\n""".format(10))
    fgrid.write("""#waveform: ricker 1.0 {} Mydipole\n""".format(paramGPRMAX.wave_freq))

    w=np.float64(max(X[0,:]))+dx/2 #largeur
    fgrid.write("""#hertzian_dipole: z {:.5f} {:.5f} {:.5f} Mydipole\n""".format(w+paramGPRMAX.d_emet, max(Y[:,0])+3/2*dy, 0.0))#TX
    fgrid.write("""#rx: {:.5f} {:.5f} {:.5f} \n""".format(w+paramGPRMAX.d_recept, max(Y[:,0]+3/2*dy), 0.0))#RX
    fgrid.write("""------------------------------------------------\n""")
    
    k=1
    VWC=np.zeros([len(Y[:,0])*len(grid_z0[0,:])+1, 1])
    for i in range(0, len(Y[:,0])) :
        for j in range(0, len(grid_z0[0,:])) :
            VWC[k-1]=grid_z0[i,j]
            k = k+1
            
    #SI EAU DOUCE INJECTEE      
    for i in materiaux:
        fgrid.write("""#material: {} {} 1.0 0.0 {}\n""".format(i, sigma, materiaux[i][0]))
        
    #SI EAU SALEE INJECTEE
    #TODO: Mettre un test pour demander à l'utilisateur quel cas il veut faire...
    #for i in materiaux:
    #    fgrid.write("""#material: {} {} 1.0 0.0 {}\n""".format(i, round(materiaux[i][1],2), materiaux[i][0]))
        
        
    #SI MATERIAUX EAU ET/OU PVC PRESENTS DANS LE MODELE
    #fgrid.write("""#material: {} {} 1.0 0.0 eau\n""".format(paramGPRMAX.eps_w,sigma)) #Donne le epsilon de l'eau
    #fgrid.write("""#material: {} {} 1.0 0.0 pvc\n""".format(paramGPRMAX.eps_pvc,sigma)) #Donne le epsilon du pvc
    fgrid.write("""------------------------------------------------\n""")
    
    #Creation du modele en entier par symetrie des resultats SWMS
    print(trace_number+1)

    d=dy/2  #hauteur
    k=1 #compteur
    count=0      
    A=np.empty([len(X[:,0])*len(X[0,:])*2, 3],dtype="float64") #Taille de la matrice = nombre de cellules du nouveau maillage*2
    A_tab={}
      
    for ii in range(0,len(Y[:,0])) :
          for jj in range(0, len(grid_z0[0,:])) :
              fgrid.write("""#box: {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {}\n""".format(w+X[ii,jj]-dx/2, d+Y[ii,jj]-dy/2, 0.0, w+X[ii,jj]+dx/2,d+Y[ii,jj]+dy/2, dl, materiaux[grid_z0[ii,jj]][0]))
              A[count,0]=w+X[ii,jj]-dx/2
              A[count,1]=d+Y[ii,jj]-dy/2
              A[count,2]=grid_z0[ii,jj]
      
              count=count+1
              fgrid.write("""#box: {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {}\n""".format(w-X[ii,jj]-dx/2, d+Y[ii,jj]-dy/2, 0.0, w-X[ii,jj]+dx/2,d+Y[ii,jj]+dy/2, dl, materiaux[grid_z0[ii,jj]][0]))
              A[count,0]=w-X[ii,jj]-dx/2
              A[count,1]=d+Y[ii,jj]-dy/2
              A[count,2]=grid_z0[ii,jj]

              k=k+1
              count=count+1
            
    #On bouche le tuyau et on interpole la fermeture du bulbe 
    fgrid.write("""------------------------------------------------\n""")
    o=np.where((A[:,0]<=w+dx/2+radius) & (A[:,0]>=w-dx/2-radius) & (A[:,1]>=etrou))
    ii=np.unique(A[o][:,1])
    blou=A[o]
    #geometry.tol=10**(-4)
    np.savetxt('bordel.csv',A,fmt='%.6f',delimiter=',')
    for i in ii:
        a=np.where(comparaison_array_number(A[o,1],i,geometry.tol))
        b=np.where((comparaison_array_number(A[:,0],w+dx/2+radius,geometry.tol))&(comparaison_array_number(A[:,1],i,geometry.tol)))
        
        for ii in a[1]:
            blou[ii,2]=A[b,2]
            A[o]=blou   
            fgrid.write("""#box: {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {}\n""".format(blou[ii,0], i,0.0, blou[ii,0]+dx,i+dy, dl, materiaux[blou[ii,2]][0]))
    #Si premier pas de temps, on ajoute une boite d'eau
    if trace_number==0 :
        q=np.where((A[:,0]<=w+dx/2+radius) & (A[:,0]>=w-dx/2-radius) & (A[:,1]>=etrou) & (A[:,1]<=etrou+h_eau))
        A[q,2]=grid_z0.max()
        fgrid.write("""#box: {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {}\n""".format(w-dx/2-radius, etrou, 0.0, w+dx/2+radius, etrou+h_eau, dl, materiaux[grid_z0.max()][0])) #materiaux[grid_z0.max()][0]
    
    #Ajout du tuyau pvc dans le trou      
    #fgrid.write("""#box: {} {} {} {} {} {} pvc\n""".format(round(w-dx/2-radius,2), round(etrou+d,2), 0.0, round(w+dx/2+radius,2), round(max(Y[:,0])*1.1,2), dl)) #ajout du tuyau tout le long
    #fgrid.write("""#box: {} {} {} {} {} {} pvc\n""".format(0.38, 0.505, 0.0, 0.43, 0.88, 0.6*dl)) #ajout du tuyau tout le long
    #o=np.where((A[:,0]<=w+dx/2+radius) & (A[:,0]>=w-dx/2-radius) & (A[:,1]>=etrou))
    #o=np.where((A[:,0]<=w+dx/2+radius) & (A[:,0]>=w-dx/2-radius) & (A[:,1]>=etrou))
    #A[o,2]=eps_pvc

    #Ajout d'eau au fond du trou
    #fgrid.write("""#box: {} {} {} {} {} {} eau\n""".format(w-dx/2-radius+0.007, etrou+d, 0.0, w+dx/2+radius-0.007, etrou+d+h_eau, dl))
    #q=np.where((A[:,0]<=w+dx/2+radius-0.007) & (A[:,0]>=w-dx/2-radius+0.007) & (A[:,1]<=etrou+d+h_eau) & (A[:,1]>=etrou+d))
    #A[q,2]=81

    #Ajout de l'air dans le trou
    #fgrid.write("""#box: {} {} {} {} {} {} free_space\n""".format(w-dx/2-radius+0.007, etrou, 0.0, w+dx/2+radius-0.007, max(Y[:,0])*2, dl))
    #s=np.where((A[:,0]<=w+dx/2+radius-0.007) & (A[:,0]>=w-dx/2-radius+0.007) & (A[:,1]<=max(Y[:,0])*1.1) & (A[:,1]>=etrou+d+h_eau))
    #A[s,2]=1
    
    #Juste pour checker la geometrie
    #fgrid.write("""#geometry_view: 0 0 f3 f4 f5 f6 f7 f8 f9 file1 c1\n""".format(dl, dl, dl))
    
    fgrid.write("""------------------------------------------------\n""")
    fgrid.write("""#messages: n\n""")
    #fgrid.write("""#geometry_view: 0 0 0 {} {} {} {} {} {} modele_vue{} n""".format(2*max(X[0,:])+dx, max(Y[:,0])*2, dl, dl, dl, dl,str(trace_number+1)))

    fgrid.close()
    
    return A
                     

