import numpy as np

def number_cells_with_BC_of_p(p):
    n_charge_imp=len(np.array(np.where(p[:,2]==1)).T)
    n_free_drainage=len(np.array(np.where(p[:,2]==-3)).T)
    return n_charge_imp + n_free_drainage

def ecriture_Selector_in(mesh, paramMVG, temps, p):
    """Fonction qui permet d'écrire le fichier Selector.in nécessaire pour lancer SWMS_2D.
    Contient surtout les param hydro et les options de SWMS_2D.
    """
    
    s = """*** BLOCK A: BASIC INFORMATION *****************************************
    'Heading'
    Example 1 - Column Test
    LUnit  TUnit  MUnit  BUnit     (units are obligatory for all input data)
     'cm'   'min'  '-'    '-'
    Kat (0:horizontal plane, 1:axisymmetric vertical flow, 2:vertical plane
    1
    MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)
      25    .001   0.5
    lWat    lChem   ChecF   ShortF  FluxF   AtmInF  SeepF  FreeD  DrainF
     t  f  f      t       t       f       f      f      f
    *** BLOCK B: MATERIAL INFORMATION **************************************
    NMat    NLay    hTab1   hTabN   NPar
      1      1      .001    200.     9
    thr     ths     tha     thm     Alfa    n       Ks      Kk      thk
     {tr}   {ts}    {tha}   {thm}   {Alfa}  {n}     {Ks}    {Kk}    {thk}
    *** BLOCK C: TIME INFORMATION ******************************************
    dt      dtMin   dtMax   DMul    DMul2   MPL
     .1    .01     10.     1.1     .7     {MPL}
    TPrint(1),TPrint(2),...,TPrint(MPL)                   (print-time array)
       {times}
    *** END OF INPUT FILE SELECTOR.IN************************************
    """.format(tr=paramMVG.tr,
               ts=paramMVG.ts,
               tha=paramMVG.tr,
               thm=paramMVG.ts,
               Alfa=paramMVG.alpha,
               n=paramMVG.n,
               Ks=paramMVG.Ks,
               Kk=paramMVG.Ks,
               thk=paramMVG.ts,
               MPL=len(temps),
               times=" ".join([str(i) for i in temps])
    )
    
    fselector=open("Selector.in","w")
    fselector.write(s)
    fselector.close()

def ecriture_Grid_in(mesh, p):
    """Fonction qui permet d'écrire le fichier Grid.in nécessaire pour lancer SWMS_2D.
    Fichier qui contient les Width...
    """
    
    Ncells = mesh.cellCount() #nombre de cellules du maillage
    Nnodes = mesh.nodeCount() #nombre de noeuds du maillage
    Nbc = number_cells_with_BC_of_p(p) #nombre de noeuds avec Boundary conditions

    t = np.zeros((Ncells, 3), dtype = 'float')
    for i in range(0,Ncells) :
        c=mesh.cell(i)
        for j in range(0,3) :
            l=c.node(j)
            f=l.id()
            t[i,j]=f
    
    a = (np.array(np.where(p[:,2]==1))).T #Noeuds ayant une charge imposée constante
    b = (np.array(np.where(p[:,2]==-3))).T #Noeuds ayant un drainage
    
    def width(p,a,b) :
        import math
        e = np.zeros((len(a), 3), dtype = 'float')
        f = np.zeros((len(b), 3), dtype = 'float')
        e[:,0] = a[:,0]
        f[:,0] = b[:,0]
        e[:,1] = p[a[:,0],0]
        f[:,1] = p[b[:,0],0]
        #On trie les noeuds par x croissant
        e=np.array(sorted(e, key=lambda colonne: colonne[1])) 
        f=np.array(sorted(f, key=lambda colonne: colonne[1]))
    
        for i in range(0,len(a[:,0]-1)) :
            if i==0 : #Premier noeud
                e[i,2] = ((math.pi)/3) * ( (e[i+1,1] + 2* e[i,1]) * (e[i+1,1]-e[i,1]))
            elif i==(len(a[:,0])-1) : #Dernier noeud
                e[i,2] = ((math.pi)/3) * ( (e[i-1,1] + 2* e[i,1]) * (e[i,1]-e[i-1,1]))
            else : #Autres noeuds
                e[i,2] = ((math.pi)/3) * ( (e[i-1,1] + 2* e[i,1]) * (e[i,1]-e[i-1,1]) + (e[i+1,1] + 2* e[i,1]) * (e[i+1,1]-e[i,1]))

                
        for i in range (0,len(b[:,0]-1)) :
            if i==0 : #Premier noeud
                f[i,2] = ((math.pi)/3) * ( (f[i+1,1] + 2* f[i,1]) * (f[i+1,1]-f[i,1]))
            elif i==(len(b[:,0])-1) : #Dernier noeud
                f[i,2] = ((math.pi)/3) * ( (f[i-1,1] + 2* f[i,1]) * (f[i,1]-f[i-1,1]))
            else : #Autre noeuds
                f[i,2] = ((math.pi)/3) * ( (f[i-1,1] + 2* f[i,1]) * (f[i,1]-f[i-1,1]) + (f[i+1,1] + 2* f[i,1]) * (f[i+1,1]-f[i,1]))
    
        return e, f 
    
    [e, f] = width(p,a,b)
    
    dim = [Nnodes, Ncells, 2, Nbc, 0]

    s1 = """*** BLOCK H: NODAL INFORMATION **************************************************
       NumNP     NumEl       IJ      NumBP     NObs
       {dim_list}
       n  Code    x      z          h       Conc      Q     M   B    Axz   Bxz   Dxz
    """.format(dim_list=" ".join([str(i) for i in dim]))

    s2 = """*** BLOCK I: ELEMENT INFORMATION ************************************************
       e   i   j   k   l   Angle  Aniz1  Aniz2  LayNum
    """

    s3="""*** BLOCK J: BOUNDARY GEOMETRY INFORMATION *************************************
        Node number array:
    """
    
    fgrid=open("Grid.in","w")

    fgrid.write(s1)

    for i in range(0,len(p)) :
        fgrid.write("""{} {} {} {} {} .00E+00 .00E+00 1 0 1 1 1 \n""".format(i + 1, int(p[i,2]), round(p[i,0],2), round(p[i,1],2), round(p[i,3],2)))
    
    #Attention, il faut enlever une ligne mais je ne vois pas comment faire...
    fgrid.write(s2)

    for i in range(0,len(t)) :
        fgrid.write("""{} {} {} {} {} 0 1 1 1 \n""".format(i + 1, int(t[i,0]+1), int(t[i,1]+1), int(t[i,2]+1), int(t[i,2]+1)))

    fgrid.write(s3)

    k = 1

    for i in range(0,len(e[:,0])) :
        fgrid.write("""{} \t""".format(int(e[i,0]+1)))
        k = k + 1
        if k==7 :
            fgrid.write("""\n""")
            k = 1
        
    k = 1

    fgrid.write("""\n""")
    if f.any() or f.size > 0:
        for i in range(0,len(f[:,0])) :
            fgrid.write("""{} \t""".format(int(f[i,0]+1)))
            k = k + 1
            if k==7 :
                fgrid.write("""\n""")
                k = 1

        k = 1 

        fgrid.write("""\n""")
    
    fgrid.write("""Width array:\n""")
    
    for i in range(0,len(e[:,0])) :
        fgrid.write("""{} \t""".format(round(e[i,2],3)))
        k = k + 1
        if k==7 :
            fgrid.write("""\n""")
            k = 1
    
    fgrid.write("""\n""")

    k = 1 

    if f.any() or f.size > 0:
        for i in range(0,len(f[:,0])) :
            fgrid.write("""{} \t""".format(round(f[i,2],3)))
            k = k + 1
            if k==7 :
                fgrid.write("""\n""")
                k=1

        fgrid.write("""\n""")
    fgrid.write("""Length:\n""")
    fgrid.write("""0\n""")

    fgrid.write("""*** END OF INPUT FILE GRID.IN *************************************************\n""")

    fgrid.close()




