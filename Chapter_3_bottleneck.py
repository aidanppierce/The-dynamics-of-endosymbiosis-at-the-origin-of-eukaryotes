
import random
import os
import numpy as np
import pandas as pd
import cProfile
from pathlib import Path
import time

def mutate(U,matrix):  
    """This counts the number of genes across the matrix, calculates how many mutations are needed and adds them to the mutation matrix"""
    
    # The number of mutations
    Nu = np.random.poisson(U, 1)

    #Locate Nu mutants in the endosymbiont matrix, this gives us the mutation_locations
    mutation_locations = random.choices(range(0,L), k = int(Nu))
    
    # this turns the numbers in mutation_locations into 3d locations and then mutates those sites
    for i in mutation_locations: 
        ng = int(i % Ng)
        ne = int((i-ng)/Ng % Ne)
        nh = int((i-(ne*Ng)-ng)/(Ne*Ng))
        if matrix[nh,ne,ng] == 1: #If a gene has already been deleted then skip
            continue
        else:
            matrix[nh,ne,ng] = 1 #If a gene hasn't been deleted then delete
    return matrix

def selection(matrix,s,theta,n): #matrix,Se,Se_neg,Se_plus,Se1,Sh,Sh_neg,Sh_plus,Sh1,kappa,Hg,Ng
    
    "We let the endosymbiont fitness be calculated at the same time as host fitness so they simulate replicating at the same time rather than one going first which skews the reproduction in favour of whatever divides first"
    ########################################
    # v10
    
    we_matrix = np.concatenate((matrix[:,:,indices[0]:indices[1]],matrix[:,:,indices[2]:indices[3]],matrix[:,:,indices[4]:indices[5]]), axis = 2)
    we = ((1-s)**np.sum(we_matrix, axis = 2))*((1+s)**np.sum(matrix[:,:,indices[3]:indices[4]], axis = 2))

    
    wh_matrix = np.concatenate((matrix[:,:,indices[1]:indices[2]],matrix[:,:,indices[3]:indices[4]],matrix[:,:,indices[4]:indices[5]]), axis = 2)
    wh = np.prod(1-(s*(np.sum(wh_matrix, axis = 1)/matrix.shape[1])**theta), axis = 1)*np.prod(1+(s*(np.sum(matrix[:,:,indices[2]:indices[3]], axis = 1)/matrix.shape[1])**theta), axis = 1)

    if n == 50000:
        Ne1 = 1
        Ne2 = 100
        new_matrix = np.zeros((Nh,Ne2,Ng))
        for i in range(0,Nh):
            #we_relative = we[i,:]/np.mean(we[i,:])
            #tot_we = np.sum(we_relative)
            #idx = np.random.choice(np.arange(Ne1), size=(Ne2), replace=True, p= we_relative/tot_we)
            
            max_location = np.where(we[i,:] == np.max(we[i,:]))[0]
            location = random.randrange(matrix.shape[1])

            for j in range(0,new_matrix.shape[1]):
                #new_matrix[i,j,:] = matrix[i,max_location[0],:]
                new_matrix[i,j,:] = matrix[i,location,:]
                

        matrix = new_matrix
    else:
        for i in range(0,Nh):
            #mean across all endosymbionts in a host and get relative fitness
            #we_relative = we[i,:]/np.mean(we[i,:])
            #tot_we = np.sum(we_relative)

            #idx = np.random.choice(np.arange(matrix.shape[1]), size=(matrix.shape[1]), replace=True, p= we_relative/tot_we)

            #########
            #Bottleneck
            max_location = np.where(we[i,:] == np.max(we[i,:]))[0]
            location = random.randrange(matrix.shape[1])

            for j in range(0,matrix.shape[1]):
                matrix[i,j,:] = matrix[i,location,:]



        
    wh_relative = wh/np.mean(wh)
    idx = np.random.choice(np.arange(Nh), size=(Nh), replace=True, p= wh_relative/np.sum(wh_relative))
    
    matrix[:,:,:] = matrix[idx,:,:]
    return matrix

def model_v10(matrix,U,s,theta,generations, repeat): #model_v9(matrix, u,Ng,Se,Se_neg,Se_plus,Se1,Sh,Sh_neg,Sh_plus,Sh1,kappa,Hg,generations, test, repeat):
    np.save(f'v10_POP_CHANGE_rand_bottlenecks_START_{repeat}_matrix', matrix)
    for n in range(0,generations):

        matrix = mutate(U,matrix)
        matrix = selection(matrix,s,theta,n)
    
        if n%1000 == 0:
            print(n)
            np.save(f'v10_POP_CHANGE_rand_bottlenecks_{n}_{repeat}_matrix', matrix)
    np.save(f'v10_POP_CHANGE_rand_bottlenecks_END_{repeat}_matrix', matrix)
    return matrix


Nh_list = [100]
Ne_list = [1]
theta_list = [1]
s_list = [0.05]
Ng = 120
u = 0.0005

generations = 100000


repeats = 50
#r_list = list(range(0,repeats+5,5))

#cProfile.run("model_v7(matrix, u,Ng,Se,Se1,Sh,Sh1,generations, test, repeat)", sort = "cumtime")
for Nh in Nh_list:
    for Ne in Ne_list:

        path = Path(f'/Users/macbook/Desktop/v10_0.0005_Ne_{Ne}_Nh_{Nh}')
        if path.is_dir():
            print('nice')
        else:
            os.mkdir(f'/Users/macbook/Desktop/v10_0.0005_Ne_{Ne}_Nh_{Nh}')
        path = f'/Users/macbook/Desktop/v10_0.0005_Ne_{Ne}_Nh_{Nh}'
        os.chdir(path)
        print(f'Nh:{Nh}, Ne:{Ne}')
        if Nh == 1000 & Ne == 1000:
            break
        for theta in theta_list:
            for s in s_list:
                print(f'Current parameters:\nNh: {Nh}, Ne: {Ne}, Ng: {Ng}, u: {u}, s: {s}, theta:{theta}')

                matrix = np.zeros((Nh,Ne,Ng))

                L = Nh * Ne * Ng
                U = u * L
                indices = np.arange(0, int(Ng+(Ng/6)), int(Ng/6))
                #The initial population is filled with genes according to the number of genes
                #start = time.time()

                parameters = pd.DataFrame({'Nh': Nh, 'Ne': Ne,'Ng':Ng,'u':u,'s': s, 'theta':theta}, index=[0]) #, 'kappa':kappa, 'Sge':Sge,'Sgh':Sgh
                parameters.to_csv(f'v10_POP_CHANGE_rand_bottlenecks_PARAMETERS.txt')

                files = [name for name in os.listdir('.') if os.path.isfile(name)]
                repeat = sum(f'v10_POP_CHANGE_rand_bottlenecks_START_' in i for i in files)
                while repeat < repeats:
                    print(f'repeat: {repeat}')
                    matrix = np.zeros((Nh,Ne,Ng))
                    model_v10(matrix,U,s,theta,generations, repeat)

                    files = [name for name in os.listdir('.') if os.path.isfile(name)]
                    repeat = sum(f'v10_POP_CHANGE_rand_bottlenecks_START_' in i for i in files)
                    
                