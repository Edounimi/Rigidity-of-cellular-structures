#!/usr/bin/env python
# coding: utf-8

# PROGRAMMA FINALE

# In[37]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import sys


# In[38]:


df = pd.read_csv('Unit_Cells.csv',",")
df = df.drop(index= [6129,8188,8190,8239,8904,8944,10843,12591,12989,13018,13283,13802,14239,14314],columns=['Unnamed: 0']) 


# Funzione per convertire coordinate e bonds dal file csv

# In[39]:


def string_toarray(dataframe,position1,position2,element):
    if element == "coordinates":
        array = dataframe.iloc[position1][position2]
        array = array.replace("[","")
        array = array.replace("]","")
        array = array.split(',')
        array = np.array(array)
        array = array.astype(float)
        dimensions = 3
        points = int(np.size(array)/dimensions)
        array = np.reshape(array,(points,dimensions))
    
        return array
    if element == "bonds":
        array = dataframe.iloc[position1][position2]
        array = array.replace("[","")
        array = array.replace("]","")
        array = array.split(',')
        array = np.array(array)
        array = array.astype(float)
        dimensions = 2
        points = int(np.size(array)/dimensions)
        array = np.reshape(array,(points,dimensions))
        array = array.astype(int)-1
        
        return array


# Funzione per plottare il reticolo

# In[40]:


def plotting_lattice(coordinates,bonds):
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter(coordinates[:,0],coordinates[:,1],coordinates[:,2], c='red', s=100)

    x= np.zeros((np.size(bonds[:,0]),2))
    y= np.zeros((np.size(bonds[:,0]),2))
    z= np.zeros((np.size(bonds[:,0]),2))


    for i in range(np.size(bonds[:,0])):
        x[i,0] = coordinates[bonds[i,0],0]
        x[i,1] = coordinates[bonds[i,1],0]
        y[i,0] = coordinates[bonds[i,0],1]
        y[i,1] = coordinates[bonds[i,1],1]
        z[i,0] = coordinates[bonds[i,0],2]
        z[i,1] = coordinates[bonds[i,1],2]

    for i in range(np.size(bonds[:,0])):
        ax.plot(x[i,:] , y[i,:], z[i,:], color='black')

    plt.show()


# Blocco di espansione reticolo

# In[41]:


def distance_from_center(coordinates,center):
    distances = np.zeros((len(coordinates[:,0]),3))
    for i in range(len(coordinates[:,0])):
        for k in range(3):
            distances[i][k] = coordinates[i][k]-center[k]
    return distances


# In[42]:


def expand_coordinates(shape,coordinates,center,traslation):
    distances = distance_from_center(coordinates,center)
    new_coordinate = np.zeros(3)
    new_center = np.zeros(3)
    for i in range(shape[0]+1):
        for k in range(shape[1]+1):
            for m in range(shape[2]+1):
                new_center[0] = center[0] + traslation*i
                new_center[1] = center[1] + traslation*k
                new_center[2] = center[2] + traslation*m
                
                for j in range(len(distances[:,0])):
                    new_coordinate[0] = new_center[0]+distances[j][0]
                    new_coordinate[1] = new_center[1]+distances[j][1]
                    new_coordinate[2] = new_center[2]+distances[j][2]
                    flag = False
                    for s in range(len(coordinates[:,0])):
                        if np.abs(coordinates[s,0] - new_coordinate[0])<np.exp(-22)  and np.abs(coordinates[s,1]-new_coordinate[1])<np.exp(-22) and np.abs(coordinates[s,2] - new_coordinate[2])<np.exp(-22) : flag = True
                    
                    if flag==False: coordinates = np.append(coordinates,[new_coordinate],axis=0)
                """
                if flag == True: 
                    print("esiste già")
                    print(new_coordinate)
                    print(" ")
                """
    """
    for l in range(len(coordinates[:,0])):
        count = 0
        for n in range(len(coordinates[:,0])):
            if coordinates[n,0] == coordinates[l,0] and coordinates[n,1] == coordinates[l,1]: count = count +1
        if count>1 : coordinates = np.delete(coordinates,(n),axis=0)  
    """
    return coordinates


# In[43]:


def expand_bonds(shape,coordinates,bonds,traslation):
    coordinate1 = np.zeros(3)
    coordinate2 = np.zeros(3)
    new_coordinate1 = np.zeros(3)
    new_coordinate2 = np.zeros(3)
    new_bond = np.zeros(2)
    new_bonds = bonds
    for s in range(len(bonds[:,0])):
        coordinate1[0] = coordinates[bonds[s,0],0]
        coordinate1[1] = coordinates[bonds[s,0],1]
        coordinate1[2] = coordinates[bonds[s,0],2]
        coordinate2[0] = coordinates[bonds[s,1],0]
        coordinate2[1] = coordinates[bonds[s,1],1]
        coordinate2[2] = coordinates[bonds[s,1],2]
        for i in range(shape[0]+1):
            for k in range(shape[1]+1):
                for j in range(shape[2]+1):
                    new_coordinate1[0] = coordinate1[0] + traslation*i
                    new_coordinate1[1] = coordinate1[1] + traslation*k
                    new_coordinate1[2] = coordinate1[2] + traslation*j
                    new_coordinate2[0] = coordinate2[0] + traslation*i
                    new_coordinate2[1] = coordinate2[1] + traslation*k
                    new_coordinate2[2] = coordinate2[2] + traslation*j
                    #print(new_coordinate1,new_coordinate2)
              
                
            
                    for l in range(len(coordinates[:,0])):
                        if np.abs(coordinates[l,0]-new_coordinate1[0])<np.exp(-22) and np.abs(coordinates[l,1]-new_coordinate1[1])<np.exp(-22) and np.abs(coordinates[l,2]-new_coordinate1[2])<np.exp(-22): new_bond[0]=l

                    for h in range(len(coordinates[:,0])):
                        if np.abs(coordinates[h,0]-new_coordinate2[0])<np.exp(-22) and np.abs(coordinates[h,1]-new_coordinate2[1])<np.exp(-22) and np.abs(coordinates[h,2]-new_coordinate2[2])<np.exp(-22): new_bond[1]=h
                
                    #print(new_bond)
                    flag = False
                
                    for v in range(len(new_bonds[:,0])):
                        if np.abs(new_bonds[v,0]-new_bond[0])<np.exp(-22) and np.abs(new_bonds[v,1]-new_bond[1])<np.exp(-22) : flag = True
                        if np.abs(new_bonds[v,1]-new_bond[0])<np.exp(-22) and np.abs(new_bonds[v,0]-new_bond[1])<np.exp(-22) : flag = True
                        if np.abs(new_bond[0]-new_bond[1])<np.exp(-22) : flag = True
                    
                    if flag==False: new_bonds = np.append(new_bonds,[new_bond],axis=0)
                
    new_bonds=new_bonds.astype(int)           
    return new_bonds


# In[44]:


def expand_lattice(coordinates,bonds,shape,traslation,center):
    new_coordinates = expand_coordinates(shape,coordinates,center,traslation)
    new_bonds = expand_bonds(shape,new_coordinates,bonds,traslation)
    return new_coordinates,new_bonds


# Periodicità dei bonds

# In[224]:


def periodic_cell(coordinates,bonds):
    
    new_coordinates = coordinates
    periodic_bonds = bonds
    nonzeropoints=np.zeros(1)
    bonds_trasformationx=np.zeros((1,2))
    bonds_trasformationy=np.zeros((1,2))
    bonds_trasformationz=np.zeros((1,2))
    
    
    #Seleziono le facce iniziali ossia quelle con coordinata nulla
    
    for i in range(len(coordinates[:,0])):
        if coordinates[i,0]!=0 and coordinates[i,1]!=0 and coordinates[i,2]!=0: nonzeropoints = np.append(nonzeropoints,i)
    
    nonzeropoints = np.delete(nonzeropoints,0)
    nonzeropoints = nonzeropoints.astype(int)
    nonzeropoints = list(set(nonzeropoints))
    new_coordinates = np.delete(new_coordinates,nonzeropoints,axis=0)
    
    #Trovo l'omologo sulle facce nei tre casi
    
    for i in range(len(new_coordinates[:,0])):
        for k in range(len(coordinates[:,0])):
            if coordinates[k,0] == new_coordinates[i,0] and coordinates[k,1] == new_coordinates[i,1] and coordinates[k,2] == new_coordinates[i,2]: point_number0 = k
        
        
        #CASO 3 COORDINATE NULLE
        if new_coordinates[i,0] <np.exp(-22) and new_coordinates[i,1] <np.exp(-22) and new_coordinates[i,2] <np.exp(-22):
            coordinate_omologuex =[2,0,0]
            coordinate_omologuey =[0,2,0]
            coordinate_omologuez =[0,0,2]
            
            Flagx = False
            Flagy = False
            Flagz = False
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuex[0] and coordinates[j,1] == coordinate_omologuex[1] and coordinates[j,2] == coordinate_omologuex[2]:
                    point_numberx = j
                    Flagx = True
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuey[0] and coordinates[j,1] == coordinate_omologuey[1] and coordinates[j,2] == coordinate_omologuey[2]:
                    point_numbery = j
                    Flagy = True
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuez[0] and coordinates[j,1] == coordinate_omologuez[1] and coordinates[j,2] == coordinate_omologuez[2]:
                    point_numberz = j
                    Flagz = True
            
            if Flagx == True:
                bonds_trasformationx = np.append(bonds_trasformationx,[[point_number0,point_numberx]],axis=0)
            
            if Flagy == True:
                bonds_trasformationy = np.append(bonds_trasformationy,[[point_number0,point_numbery]],axis=0)
            
            if Flagz == True:
                bonds_trasformationz = np.append(bonds_trasformationz,[[point_number0,point_numberz]],axis=0)
            
            
        #CASI 2 COORDINATE NULLE
        
        elif new_coordinates[i,0] <np.exp(-22) and new_coordinates[i,1] <np.exp(-22):
            coordinate_omologuex =[2,0,new_coordinates[i,2]]
            coordinate_omologuey =[0,2,new_coordinates[i,2]]
            
            Flagx = False
            Flagy = False
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuex[0] and coordinates[j,1] == coordinate_omologuex[1] and coordinates[j,2] == coordinate_omologuex[2]:
                    point_numberx = j
                    Flagx = True
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuey[0] and coordinates[j,1] == coordinate_omologuey[1] and coordinates[j,2] == coordinate_omologuey[2]:
                    point_numbery = j
                    Flagy = True
            
            if Flagx == True:
                bonds_trasformationx = np.append(bonds_trasformationx,[[point_number0,point_numberx]],axis=0)
            
            if Flagy == True:
                bonds_trasformationy = np.append(bonds_trasformationy,[[point_number0,point_numbery]],axis=0)
        
        
        elif new_coordinates[i,1] <np.exp(-22) and new_coordinates[i,2] <np.exp(-22):
            coordinate_omologuez =[new_coordinates[i,0],0,2]
            coordinate_omologuey =[new_coordinates[i,0],2,0]
            
            Flagz = False
            Flagy = False
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuez[0] and coordinates[j,1] == coordinate_omologuez[1] and coordinates[j,2] == coordinate_omologuez[2]:
                    point_numberz = j
                    Flagz = True
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuey[0] and coordinates[j,1] == coordinate_omologuey[1] and coordinates[j,2] == coordinate_omologuey[2]:
                    point_numbery = j
                    Flagy = True
            
            if Flagz == True:
                bonds_trasformationz = np.append(bonds_trasformationz,[[point_number0,point_numberz]],axis=0)
            
            if Flagy == True:
                bonds_trasformationy = np.append(bonds_trasformationy,[[point_number0,point_numbery]],axis=0)
            
            
        elif new_coordinates[i,0] <np.exp(-22) and new_coordinates[i,2] <np.exp(-22):
            coordinate_omologuex =[2,new_coordinates[i,1],0]
            coordinate_omologuez =[0,new_coordinates[i,1],2]
            
            Flagx = False
            Flagz = False
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuex[0] and coordinates[j,1] == coordinate_omologuex[1] and coordinates[j,2] == coordinate_omologuex[2]:
                    point_numberx = j
                    Flagx = True
            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuez[0] and coordinates[j,1] == coordinate_omologuez[1] and coordinates[j,2] == coordinate_omologuez[2]:
                    point_numberz = j
                    Flagz = True
            
            if Flagx == True:
                bonds_trasformationx = np.append(bonds_trasformationx,[[point_number0,point_numberx]],axis=0)
                
            if Flagz == True:
                bonds_trasformationz = np.append(bonds_trasformationz,[[point_number0,point_numberz]],axis=0)
            
        #CASI AD UNA COORDINATA NULLA
        
        elif new_coordinates[i,0] <np.exp(-22):
            coordinate_omologuex =[2,new_coordinates[i,1],new_coordinates[i,2]]

            Flagx = False

            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuex[0] and coordinates[j,1] == coordinate_omologuex[1] and coordinates[j,2] == coordinate_omologuex[2]:
                    point_numberx = j
                    Flagx = True
            
            if Flagx == True:
                bonds_trasformationx = np.append(bonds_trasformationx,[[point_number0,point_numberx]],axis=0)
        
        elif new_coordinates[i,1] <np.exp(-22):
            coordinate_omologuey =[new_coordinates[i,0],2,new_coordinates[i,2]]

            Flagy = False

            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuey[0] and coordinates[j,1] == coordinate_omologuey[1] and coordinates[j,2] == coordinate_omologuey[2]:
                    point_numbery = j
                    Flagy = True
            
            if Flagy == True:
                bonds_trasformationy = np.append(bonds_trasformationy,[[point_number0,point_numbery]],axis=0)
        
        
        elif new_coordinates[i,2] <np.exp(-22):
            coordinate_omologuez =[new_coordinates[i,0],new_coordinates[i,1],2]

            Flagz = False

            
            for j in range(len(coordinates[:,0])):
                if coordinates[j,0] == coordinate_omologuez[0] and coordinates[j,1] == coordinate_omologuez[1] and coordinates[j,2] == coordinate_omologuez[2]:
                    point_numberz = j
                    Flagz = True
            
            if Flagz == True:
                bonds_trasformationz = np.append(bonds_trasformationz,[[point_number0,point_numberz]],axis=0)
    """
    print("Bonds prima")            
    print(bonds)
    """
    #INIZIO LE TRASFORMAZIONI

    bonds_trasformationx = np.delete(bonds_trasformationx,0,axis=0)
    bonds_trasformationy = np.delete(bonds_trasformationy,0,axis=0)
    bonds_trasformationz = np.delete(bonds_trasformationz,0,axis=0)
    """
    print(bonds_trasformationx)
    print(bonds_trasformationy)
    print(bonds_trasformationz)
    """
    bonds_real = np.zeros((len(bonds[:,0]),2))
    for i in range(len(bonds[:,0])):
        bonds_real[i,1] = bonds[i,1]
        bonds_real[i,0] = bonds[i,0]
    
    #TRASFORMAZIONE LUNGO X
    
    for i in range(len(bonds_trasformationx[:,0])):
        for s in range(len(periodic_bonds[:,0])):
            if np.abs(periodic_bonds[s,0]-bonds_trasformationx[i,1])<np.exp(-22):
                periodic_bonds[s,0]= bonds_trasformationx[i,0]
                
            elif np.abs(periodic_bonds[s,1]-bonds_trasformationx[i,1])<np.exp(-22):
                periodic_bonds[s,1]= bonds_trasformationx[i,0]
    
    #TRASFORMAZIONE LUNGO Y
    
    for i in range(len(bonds_trasformationy[:,0])):
        for s in range(len(periodic_bonds[:,0])):
            if np.abs(periodic_bonds[s,0]-bonds_trasformationy[i,1])<np.exp(-22):
                periodic_bonds[s,0]= bonds_trasformationy[i,0]
                
            elif np.abs(periodic_bonds[s,1]-bonds_trasformationy[i,1])<np.exp(-22):
                periodic_bonds[s,1]= bonds_trasformationy[i,0]
    
    #TRASFORMAZIONE LUNGO Z
    
    for i in range(len(bonds_trasformationz[:,0])):
        for s in range(len(periodic_bonds[:,0])):
            if np.abs(periodic_bonds[s,0]-bonds_trasformationz[i,1])<np.exp(-22):
                periodic_bonds[s,0]= bonds_trasformationz[i,0]
                
            elif np.abs(periodic_bonds[s,1]-bonds_trasformationz[i,1])<np.exp(-22):
                periodic_bonds[s,1]= bonds_trasformationz[i,0]
    
    """
    print("bonds real")
    print(bonds_real)
    print(periodic_bonds)
    
    """
    #TOLGO I BONDS NON COINVOLTI NELLA TRASFORMAZIONE
    
    
    new_periodic_bonds_indices = np.array(0)
    
    for i in range(len(periodic_bonds[:,0])):
        if np.abs(periodic_bonds[i,0]-bonds_real[i,0])<np.exp(-22) and np.abs(periodic_bonds[i,1]-bonds_real[i,1])<np.exp(-22): 
            new_periodic_bonds_indices = np.append(new_periodic_bonds_indices,i)
            
    new_periodic_bonds_indices = np.delete(new_periodic_bonds_indices,0)
    new_periodic_bonds_indices = list(set(new_periodic_bonds_indices))
    periodic_bonds = np.delete(periodic_bonds,new_periodic_bonds_indices,axis=0)
   

    
    #VERIFICHIAMO CHE NON CI SIANO BONDS RIPETUTI
    
    periodic_bonds_temp = np.array(0)
    for i in range(int(len(periodic_bonds[:,0])-1)):
        for m in range(i+1,len(periodic_bonds[:,0])):
            if periodic_bonds[m,0] == periodic_bonds[i,0] and periodic_bonds[m,1] == periodic_bonds[i,1] or periodic_bonds[m,0] == periodic_bonds[i,1] and periodic_bonds[m,1] == periodic_bonds[i,0]: periodic_bonds_temp = np.append(periodic_bonds_temp,m)
    
    
    periodic_bonds_temp = np.delete(periodic_bonds_temp,0)
    periodic_bonds_temp = list(set(periodic_bonds_temp))
    periodic_bonds = np.delete(periodic_bonds,periodic_bonds_temp,axis=0)
    
    
    return periodic_bonds, bonds_real.astype(int)
    


# Calcolo delle matrici di adiacenza

# In[225]:


def adjacency_matrix(coordinates,bonds,periodic_bonds=None,PBC = False):
    N= len(coordinates[:,0])
    adjacency_matrix = np.zeros((N,N))
    
    for i in range(len(bonds[:,0])):
        r = bonds[i,0]
        c = bonds[i,1]
        
        adjacency_matrix[r,c]=1
        adjacency_matrix[c,r]=1
    
    if PBC == True:
        for i in range(len(periodic_bonds[:,0])):
            r_p = periodic_bonds[i,0]
            c_p = periodic_bonds[i,1]
            adjacency_matrix[r_p,c_p]+=1
            adjacency_matrix[c_p,r_p]+=1
            
    return adjacency_matrix


# Calcolo della matrice hessiana

# In[226]:


def Hessian_super_element_non_diagonal(AM,coordinates,i,j):
    H_nd = np.zeros((3,3))
    
    distance = np.sqrt(np.power(coordinates[i,0]-coordinates[j,0],2)+np.power(coordinates[i,1]-coordinates[j,1],2)+np.power(coordinates[i,2]-coordinates[j,2],2))
    for k in range(3):
        for l in range(3):
            H_nd[k,l] = -AM[i,j]*(coordinates[i,k]-coordinates[j,k])*(coordinates[i,l]-coordinates[j,l])/np.power(distance,2)
    
    #print(H_nd)
    return H_nd

def Hessian_super_element_diagonal(AM,coordinates,i):
    H_d = np.zeros((3,3))
    
    supersum = 0 
    
    for j in range(AM.shape[0]):
        if j!=i:
            supersum = supersum + Hessian_super_element_non_diagonal(AM,coordinates,i,j)
    
    H_d = -supersum
    #print(H_d)
    return H_d
    
def Hessian_Matrix_Theory(AM,coordinates):
    
    H = np.zeros((3*AM.shape[0],3*AM.shape[0]))
    
    for i in range(AM.shape[0]):
        for j in range(AM.shape[0]):
            if j!=i:
                H[i*3:i*3+3,j*3:j*3+3]=Hessian_super_element_non_diagonal(AM,coordinates,i,j)
            if j == i:
                H[i*3:i*3+3,j*3:j*3+3]=Hessian_super_element_diagonal(AM,coordinates,j)
    return H


# In[227]:


def count_zeroeigenvalues(eigenvalues):
    k = 0
    for i in range(np.size(eigenvalues)):
        if np.abs(eigenvalues[i].real) < np.exp(-22) : 
            k = k + 1
    
    return k-6


# ANALISI COMPLETA

# In[244]:


def lattice_analysis_complete(dataframe,cell_number,base_cell = False,plot=False,PBC = False):
    
    
    #names = dataframe["Name"].tolist()
    
    coordinates = string_toarray(dataframe,cell_number,27,element="coordinates")
    bonds = string_toarray(dataframe,cell_number,28,element="bonds")
    
    """
    print("Analisi della cella numero " + str(cell_number))
    print("Nome cella: " + names[cell_number])
    """
    if base_cell == True:
        print("Cella di base")
        plotting_lattice(coordinates,bonds)

    shape=np.array([1,1,1])
    center = np.array([0.5,0.5,0.5])
    traslation = 1
    
    #print("Espando il reticolo..")
    coordinates_expanded, bonds_expanded = expand_lattice(coordinates,bonds,shape,traslation,center)
    AM = adjacency_matrix(coordinates_expanded,bonds_expanded,PBC = False)
    H = Hessian_Matrix_Theory(AM,coordinates_expanded)
    Eigenvalues, Eigenvectors = np.linalg.eig(H)
    
    """
    print("Matrice di adiacenza")
    print(AM)
    print("Matrice hessiana")
    print(H)
    print("Autovalori del reticolo espanso: ")
    print(Eigenvalues.real)
    """
    
    NT_zeroEigenvalues = count_zeroeigenvalues(Eigenvalues)
    
    #print("Autovalori non triviali nulli : " + str(NT_zeroEigenvalues))
    
    if plot == True:
        print("Rappresentazione grafica del reticolo")
        plotting_lattice(coordinates_expanded,bonds_expanded)
        
    
    if PBC == True:
        #print("Applico le PBC..")
        periodic_bonds,bonds_real = periodic_cell(coordinates_expanded,bonds_expanded)
        AM_P = adjacency_matrix(coordinates_expanded,bonds_real,periodic_bonds,PBC = True)
        H_P = Hessian_Matrix_Theory(AM_P,coordinates_expanded)
        Eigenvalues_P, Eigenvectors_P = np.linalg.eig(H_P)
        
        """
        print("Matrice di adiacenza con PBC")
        print(AM_P)
        print("Matrice hessiana con PBC")
        print(H_P)
        print("Autovalori del reticolo espanso con PBC: ")
        print(Eigenvalues_P.real)
        """
        NT_zeroEigenvalues_P = count_zeroeigenvalues(Eigenvalues_P)
        #print("Autovalori non triviali nulli con PBC : " + str(NT_zeroEigenvalues_P))
        #print("   ")
        Number_Eigenvalues = len(Eigenvalues_P)
        return Number_Eigenvalues,NT_zeroEigenvalues, NT_zeroEigenvalues_P
        
    if PBC == False:
        print("   ")
        return NT_zeroEigenvalues



# In[236]:


i=int(sys.argv[1])
Eigendata=lattice_analysis_complete(df,i,base_cell = False,PBC = True,plot =False)
print(Eigendata)

data= {'Cell number':[i],'Eigenvalues number':[Eigendata[0]],
           'Non trivial zero Eigenvalues': [Eigendata[1]],
           'Non trivial zero Eigenvalues con PBC': [Eigendata[2]]}

df_final = pd.DataFrame(data, columns=['Cell number','Eigenvalues number','Non trivial zero Eigenvalues', 'Non trivial zero Eigenvalues con PBC'])

title = "cellnumber_"+str(i)+".csv"
csv_data = df_final.to_csv(title)



