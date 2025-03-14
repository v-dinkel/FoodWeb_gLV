# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 17:05:34 2022

@author: vdinkel
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from scipy import integrate
import networkx as nx
import networkx.algorithms.community as nx_comm

def msum(A):
    return int(np.sum(np.sum(A)))

def AND(A, B):
    return B.where(A == 1, 0)

def fillDiagonal(dataFrame, value = 0):
    for index in dataFrame.columns:
        dataFrame.at[index,index] = value
    return dataFrame

def binNorm(matrix):
    retM = matrix.where(matrix == 0, 1)
    return retM

def ESABO(x1, x2):
    # ESABO gives H and z score of the logical and from x1 and x2
    xij = np.logical_and(x1,x2).astype(int)
    p1 = np.sum(xij) / len(xij)
    
    p1a = np.sum(x1) / len(x1)
    p1b = np.sum(x2) / len(x2)
    
    std = np.sqrt(  (p1a * p1b * (1-(p1a*p1b))) / len(x1)  )
    mean = p1a * p1b

    z_score = (p1 - mean) / std
    
    H_real = mean
    return H_real, z_score

def corrESABOM(matrix, thresh):
    # computes ESABO network from given abundance matrix with thresh being the threshold for the z-score.
    ESABO_A = np.zeros((matrix.shape[1],matrix.shape[1]))
    for i in range(0,matrix.shape[1]):
        for j in range(0,matrix.shape[1]):
            if i==j:
                ESABO_A[i][j] = 0
            else:
                arr1 = matrix[matrix.columns[i]].values # abundance vector of taxon i
                arr2 = matrix[matrix.columns[j]].values # abundance vector of taxon j
                if (0 in arr1 and 1 in arr1 and 0 in arr2 and 1 in arr2):
                    H_real, z_score = ESABO(arr1, arr2)
                    if np.isnan(z_score):
                        import pdb; pdb.set_trace()
                    ESABO_A[i][j] = z_score
                else:
                    ESABO_A[i][j] = 0
                
    ESABO_A = pd.DataFrame(ESABO_A)
    ESABO_A.columns = matrix.columns
    ESABO_A.index = matrix.columns
    
    ESABO_A = ESABO_A.where(abs(ESABO_A) > thresh, 0) # z score must be > threshold (set to 1.3)

    return ESABO_A

def dX_dt_gLV(M, t=0):
    # glV funktion
    return r*M * (1 - np.dot(A, M) / k)

def getX0(N):
    # initial definition of parameters for lotka volterra simulation.
    
    # M0 => intial population counts
    # initial abundances of the species were drawn from a uniform distribution ranging between 10 and 100.
    x0 = []
    for i in range(0,N):
        x0.append(random.randint(10,100))
    x0 = np.array(x0)
    return x0

def getK(N):
    # CARRYING CAPACITY beta distribution, but also normal distribution.
    # scaled to range from 1 to 100
    # returns a vector of carrying capacities
    k = []
    for i in range(0,N):
        k.append(random.randint(1,100))
    k = np.array(k)
    return k

def getR(N):
    # B => growth rates
    # Growth rates were assigned to each species from a uniform distribution from 0 to 1 so that all species were capable of positive growth
    # returns a vector of growth rates
    r = []
    for i in range(0,N):
        r.append(random.randint(1,100)/100)
    r = np.array(r)
    return r

def modularity(M):
    if type(M) == type(pd.DataFrame()):
        G = nx.from_numpy_matrix(np.matrix(M))
    else:
        G = M
    comms = nx_comm.greedy_modularity_communities(G, resolution=1.0)
    return float(round(nx_comm.modularity(G, comms),2)), comms

def step1_createLoadA(workdir, base_graph_path):
    # this function creates a weighted interaction matrix
    # read simulation network topology
    sim_A = pd.read_csv(base_graph_path, delimiter=",", header=0, index_col=0)
    sim_A.index = sim_A.columns
    
    # create matrix with nxn random values between 0.5 and 1.0 as interactions strengths
    N = sim_A.shape[0]
    multM = np.random.uniform(low=.5, high=1.0, size=(N,N))
    
    # get random nxn matrix of [-1,1] which encodes the interaction types. Symmetrical means no predation etc.
    randInts = np.asmatrix(np.reshape(np.random.choice([-1,1], N*N), (N,N) ))
    randInts_symm = np.triu(randInts, k=1) + np.triu(randInts, k=1).T
    
    # this multiplication joins the matrices of topology, interaction strength and interaction type.
    new_sim_A = sim_A * multM * randInts
    new_sim_A = new_sim_A.replace([-0], 0)
    new_sim_A = round(new_sim_A,1)
    
    new_sim_A.index = pd.Index([k for k in range(1,N+1)]) 
    new_sim_A.columns = pd.Index([k for k in range(1,N+1)]) 
    # set the diagonal to 1
    new_sim_A = fillDiagonal(new_sim_A, value = 1.0)
    
    new_sim_A.index = new_sim_A.columns
    return new_sim_A, sim_A

def filtMatrices(abunds, new_sim_A, base_sim_A):
    # this function identifies species with abundance sum = 0. Some inference methods have this requirement.
    abunds.columns = base_sim_A.columns
    new_sim_A.columns = base_sim_A.columns
    new_sim_A.index = base_sim_A.index

    badcols = []
    badindeces = []
    i = 0

    for col in abunds.columns:
        if sum(abunds[col]) == 0:
            badcols.append(col)
            badindeces.append(i)
        i += 1
    filt_abunds = abunds.drop(abunds.columns[badindeces],axis = 1)
    
    filt_new_sim_A = new_sim_A.drop(new_sim_A.columns[badindeces],axis = 1)
    filt_new_sim_A = filt_new_sim_A.drop(filt_new_sim_A.index[badindeces],axis = 0)
    
    filt_base_sim_A = base_sim_A.drop(base_sim_A.columns[badindeces],axis = 1)
    filt_base_sim_A = filt_base_sim_A.drop(filt_base_sim_A.index[badindeces],axis = 0)
    
    return filt_abunds, filt_base_sim_A, filt_new_sim_A, badcols


# read the snakemake config variables and input network
workdir = snakemake.config['workdir']
nsimulations = snakemake.config['nsimulations']
base_graph_dir = snakemake.input[0]

# amount of samples (attractor states) per simulation
n_samples = 100
t = np.linspace(0.0, 100.0,  100)

nrun = 0
while nsimulations >0:
    print ("--- nsimulations: ",nsimulations)
    new_sim_A, base_sim_A = step1_createLoadA(workdir, base_graph_dir)
    N = new_sim_A.shape[0]
    
    # perform single gLV simulation
    A = new_sim_A
    x0 = getX0(N)
    k = getK(N)
    r = getR(N)
    X, infodict = integrate.odeint(dX_dt_gLV, x0, t, full_output=True)
    
    testV = np.matrix(abs(X)).flatten()
    testV = np.array(testV)[0]
    
    if max(testV) > 10000 or sum(X[99]) == 0:
        # the highest abundance value must not exceed 10000 (indicates an artefact as exponential peak in simulation)
        # the sum of the attractors must not be 0
        # failed simulation indicates unstable state of interaction matrix, therefore it needs to be rebuild
        print ("-- FAILED SIMULATION")
    else:
        print ("-- PASSED SIMULATION")
        # simulation matrix has passed, so now the attractor vector can be created
        # create and save the simulation plot
        plt.plot(X)
        plt.title("Single Simulation")
        plt.xlabel("Time")
        plt.ylabel("Abundance")
        plt.savefig(snakemake.output[nrun].replace(".csv","_sim_.jpg"))
        plt.clf()
        
        # perform defined amount of simulations, store last values (attractor) to abundance matrix
        abunds = []
        for j in range(0, n_samples):
            k = getK(N)
            r = getR(N)
            X, infodict = integrate.odeint(dX_dt_gLV, x0, t, full_output=1)
            abunds.append(X[99]) #take only the last values (attractors)
        
        if (max(np.array(abs(np.matrix(abunds).flatten()))[0]) > 10000) or sum(pd.DataFrame(abunds).T.sum()>1) <20:
            # one more check of the final attractor values and the matrix to be within the range of >20 and <10000.
            print ("--- FAILED ATTRACTORS") 
        else:
            print ("--- PASSED ATTRACTORS")
            # plot and save attractor matrix
            plt.plot(abunds)
            plt.title("Simulated Attractor Abundances")
            plt.xlabel("Simulation")
            plt.ylabel("Attracted Abundance")
            plt.savefig(snakemake.output[nrun].replace(".csv","_abunds_.jpg"))
            plt.clf()
    
            # abundances are only kept if they are > 1 (no 0.x abundances)
            abunds = pd.DataFrame(abunds)
            abunds = abunds.where(abunds >= 1, 0)
    
            # remove zero sums
            filt_abunds, filt_base_A, filt_new_sim_A, badcols = filtMatrices(abunds, fillDiagonal(new_sim_A,0) , base_sim_A) #fillDiagonal(new_sim_A, 0)
            print ("-- REMOVED ", str(len(badcols)), " Indeces")
    
            # save this simulation files
            pd.DataFrame.to_csv(abunds, snakemake.output[nrun], index=True, sep=",")
            
            pd.DataFrame.to_csv(pd.DataFrame(new_sim_A), snakemake.output[nrun].replace(".csv","_new_sim_A.csv"))
            pd.DataFrame.to_csv(pd.DataFrame(base_sim_A), snakemake.output[nrun].replace(".csv","_base_sim_A.csv"))
            
            pd.DataFrame.to_csv(pd.DataFrame(filt_abunds), snakemake.output[nrun].replace(".csv","_filt_abunds.csv"))
            pd.DataFrame.to_csv(pd.DataFrame(filt_base_A), snakemake.output[nrun].replace(".csv","_filt_base_A.csv"))
            pd.DataFrame.to_csv(pd.DataFrame(filt_new_sim_A), snakemake.output[nrun].replace(".csv","_filt_new_sim_A.csv"))
            
            #infer ESABO network
            matrix = binNorm(filt_abunds)
            esaboM = corrESABOM(matrix, 1.3) #z-score threshold is set to 1.3.
            absEsabo = abs(esaboM.round().astype(int))
            absEsabo = absEsabo.where(absEsabo == 0, 1)
            
            #save ESABO network
            pd.DataFrame.to_csv(absEsabo, snakemake.output[nrun].replace("/abundances/","/networks/").replace(".csv","_esabo.csv"))
            
            # only when simulations were successful, continue counting (backwards) until the nsimulation amount is reached
            nsimulations -= 1
            nrun += 1