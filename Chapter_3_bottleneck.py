
import random
import os
import numpy as np
import pandas as pd
import cProfile
from pathlib import Path
import time

def mutate(matrix):  
    """This counts the number of genes across the matrix, calculates how many mutations are needed and adds them to the mutation matrix"""
    
    # The number of mutations
    L = Nh * matrix.shape[1] * Ng #Define the total number of loci
    U = u * L #Define the total number of mutations per generation (on average)
    Nu = np.random.poisson(U, 1)

    #Locate Nu mutants in the endosymbiont matrix, this gives us the mutation_locations
    mutation_locations = random.choices(range(0,L), k = int(Nu))
    
    # this turns the numbers in mutation_locations into 3d locations and then mutates those sites
    for i in mutation_locations: 
        ng = int(i % Ng)
        ne = int((i-ng)/Ng % matrix.shape[1])
        nh = int((i-(ne*Ng)-ng)/(matrix.shape[1]*Ng))
        if matrix[nh,ne,ng] == 1: #If a gene has already been deleted then skip
            continue
        else:
            matrix[nh,ne,ng] = 1 #If a gene hasn't been deleted then delete
    return matrix

def selection(matrix,s,Ne,n): #matrix,Se,Se_neg,Se_plus,Se1,Sh,Sh_neg,Sh_plus,Sh1,kappa,Hg,Ng
    
    "We let the endosymbiont fitness be calculated at the same time as host fitness so they simulate replicating at the same time rather than one going first which skews the reproduction in favour of whatever divides first"
    ########################################
    
    we_matrix = np.concatenate((matrix[:,:,indices[0]:indices[1]],matrix[:,:,indices[2]:indices[3]],matrix[:,:,indices[4]:indices[5]]), axis = 2)
    we = ((1-s)**np.sum(we_matrix, axis = 2))*((1+s)**np.sum(matrix[:,:,indices[3]:indices[4]], axis = 2))

    
    wh_matrix = np.concatenate((matrix[:,:,indices[1]:indices[2]],matrix[:,:,indices[3]:indices[4]],matrix[:,:,indices[4]:indices[5]]), axis = 2)
    wh = np.prod(1-(s*(np.sum(wh_matrix, axis = 1)/matrix.shape[1])**theta), axis = 1)*np.prod(1+(s*(np.sum(matrix[:,:,indices[2]:indices[3]], axis = 1)/matrix.shape[1])**theta), axis = 1)

    if n == generations/2:
        Ne1 = Ne
        Ne2 = 100
        new_matrix = np.zeros((Nh,Ne2,Ng))
        for i in range(0,Nh):
            
            #if bottleneck == 'None':
            we_relative = we[i,:]/np.mean(we[i,:])
            tot_we = np.sum(we_relative)
            location = np.random.choice(np.arange(Ne1), size=(Ne2), replace=True, p= we_relative/tot_we)
                
            #if bottleneck == 'Fittest':
                ##only the fittest endosymbiont reproduces
            #fittest = np.where(we[i,:] == np.max(we[i,:]))[0]
            #max_location = random.randrange(len(fittest))
            #if bottleneck == 'Random':
                ##Completely random reproduction
            #rand_location = random.randrange(matrix.shape[1])

            #for j in range(0,new_matrix.shape[1]):
            #For the fittest endosymbiont
            #new_matrix[i,j,:] = matrix[i,max_location[0],:]

            new_matrix[i,:,:] = matrix[i,location,:]
                

        matrix = new_matrix
    else:
        for i in range(0,Nh):
            #mean across all endosymbionts in a host and get relative fitness
            we_relative = we[i,:]/np.mean(we[i,:])
            tot_we = np.sum(we_relative)

            #########
            #Bottleneck reproduction
            #normal likelihood based replication
            location = np.random.choice(np.arange(matrix.shape[1]), size=(matrix.shape[1]), replace=True, p= we_relative/tot_we)

            #only the fittest endosymbiont reproduces
            #fittest = np.where(we[i,:] == np.max(we[i,:]))[0]
            #location = random.randrange(len(fittest))
            
            #Completely random reproduction
            #location = random.randrange(matrix.shape[1])

            #for j in range(0,matrix.shape[1]):
            matrix[i,:,:] = matrix[i,location,:]
        
    wh_relative = wh/np.mean(wh)
    idx = np.random.choice(np.arange(Nh), size=(Nh), replace=True, p= wh_relative/np.sum(wh_relative))
    
    matrix[:,:,:] = matrix[idx,:,:]
    return matrix

def model_bottleneck(matrix,s,theta,generations, repeat): 
    np.save(f'{bottleneck_type}_bottleneck_Ne{Ne}_Nh{Nh}_u{u}_s{s}_theta{theta}_{Ng}_START_{repeat}_matrix', matrix)
    for n in range(0,generations):

        matrix = mutate(matrix)
        matrix = selection(matrix,s,Ne,n)
    
        if n%1000 == 0:
            print(n)
            np.save(f'{bottleneck_type}_bottleneck_Ne{Ne}_Nh{Nh}_u{u}_s{s}_theta{theta}_{Ng}_{n}_{repeat}_matrix', matrix)
    np.save(f'{bottleneck_type}_bottleneck_Ne{Ne}_Nh{Nh}_u{u}_s{s}_theta{theta}_{Ng}_{generations}_{repeat}_matrix', matrix)
    return matrix



#The parameters you want to test, these are as lists so they can be looped through 
Nh_list = [100] #The number of hosts
Ne_list = [1] #The number of endosymbionts
theta_list = [1] #The value of the epistasis coefficient
s_list = [0.05] #The value of the selection coefficient
Ng_list = [120] #The total number of genes in a starting endosymbiont (120: 40 per gene)
u_list = [0.0005] #The per gene rate of mutation

generations = 100000 #How many generations will the simulations go on for

#Determine how many terminals you want to run this over for parallel running of code
terminal_number = input('what terminal number is this?')
repeats = np.arange(5*(int(terminal_number)-1),5*int(terminal_number))
bottleneck_type = 'fittest'

#This is for looping through the simulations so it can be done continuously
for Nh in Nh_list:
    for Ne in Ne_list:
        for theta in theta_list:
            for Ng in Ng_list:
                #Name the directory you want the code to go to, this is an example directory
                path = Path(f'/Users/macbook/Desktop/bottleneck')
                if path.is_dir():
                    print('nice')
                else:
                    os.mkdir(f'/Users/macbook/Desktop/bottleneck')
                path = f'/Users/macbook/Desktop/bottleneck'
                os.chdir(path)
                print(f'Nh:{Nh}, Ne:{Ne}')
                if Nh == 1000 & Ne == 1000:
                    break
                for s in s_list:
                    for u in u_list:
                        print(f'Current parameters:\nNh: {Nh}, Ne: {Ne}, Ng: {Ng}, u: {u}, s: {s}, theta:{theta}')
                        
                        L = Nh * Ne * Ng #Define the total number of loci
                        U = u * L #Define the total number of mutations per generation (on average)
                        indices = np.arange(0, int(Ng+(Ng/6)), int(Ng/6))

                        for repeat in repeats:
                            print(f'repeat: {repeat}')
                            matrix = np.zeros((Nh,Ne,Ng))
                            model_bottleneck(matrix,s,theta,generations,repeat)
print('Done')
