# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 17:05:34 2022

@author: vdinkel
"""

#from numpy import *
import numpy as np
import pandas as pd
import pylab as p
import matplotlib.pyplot as plt
import random
from numpy.random import choice
from scipy import integrate
from scipy import stats
import networkx as nx
import networkx.algorithms.community as nx_comm

def msum(A):
    return int(np.sum(np.sum(A)))

def AND(A, B):
    return B.where(A == 1, 0)

def OR(A, B):
    C = abs(A) + abs(B)
    return C.replace(2, 1)

def fillDiagonal(dataFrame, value = 0):
    for index in dataFrame.columns:
        dataFrame.at[index,index] = value
    return dataFrame

def getOwnTPs(trueA, predA, threshold):
    tp = 0;fp = 0; tn = 0; fn = 0
    
    i = 0
    for predArr in predA.values:
        j = 0
        for predVal in predArr:
            predVal = abs(predVal)
            if predVal < threshold:
                predVal = 0
            if predVal > 0 and trueA.values[i][j] > 0:
                tp += 1
            if predVal > 0 and trueA.values[i][j] == 0:
                fp += 1
            if predVal == 0 and trueA.values[i][j] == 0:
                tn += 1
            if predVal == 0 and trueA.values[i][j] > 0:
                fn += 1
            j+=1
        i += 1
    if tp+fn != 0:
        tpr = tp/(tp+fn)
    else:
        tpr = 0
    if fp+tn != 0:
        fpr = fp/(fp+tn)
    else:
        fpr = 0
    if tp+fp+tn+fp != 0:
        acc = (tp+tn) / (tp+fp+tn+fp)
    else:
        acc = 0
    return {"tp": tp, "fp":fp, "tn": tn, "fn": fn, "tpr": tpr, "fpr": fpr, "acc": acc}

def compareNetworks(baseA, testB):
    intersect_all = AND(abs(baseA), abs(testB))
    try:
        JI_all = (msum(intersect_all)) / min(msum(abs(baseA)), msum(abs(testB))) #(msum(abs(baseA)) + msum(abs(testB)) - msum(abs(intersect_all)))
    except:
        JI_all = 0
    return JI_all


def binNorm(matrix):
    retM = matrix.where(matrix == 0, 1)
    return retM

def ESABO(x1, x2):
    # ESABO test
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
    # columns = taxa/features; rows = samples/observations
    
    ESABO_A = np.zeros((matrix.shape[1],matrix.shape[1])) # matrix[y][x]
    #ESABO_p = np.zeros((matrix.shape[0],matrix.shape[0])) # matrix[y][x]
    
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
                    #print ("not compatible with esabo")
                    ESABO_A[i][j] = 0
                
    ESABO_A = pd.DataFrame(ESABO_A)
    ESABO_A.columns = matrix.columns
    ESABO_A.index = matrix.columns
    
    ESABO_A = ESABO_A.where(abs(ESABO_A) > thresh, 0) # z score must be > 2
    #ESABO_A[ESABO_A != 0] = 1
    
    return ESABO_A#.astype(int)

def dX_dt_gLV(M, t=0):
    '''
    newM = []
    for i in range(0,N):
        sum = 0
        for j in range (0, N):
            sum += A[i][j] * M[j]
        newM.append(max(0, int(M[i]*B[i] + M[i]*sum)))
    return np.array(newM)
    #return array([ a*X[0] -   b*X[0]*X[1] , -c*X[1] + d*b*X[0]*X[1] ])
    '''
    
    # ri * xi * (1 - Ax/ki)
    
    #import pdb; pdb.set_trace()
    #return np.where(M*B+np.dot(A,(M.reshape(1,len(M)).T)).T[0] < 0.0,0.0, M*B+np.dot(A,(M.reshape(1,len(M0)).T)).T[0])
    #B*M+np.dot(A,(M.reshape(1,len(M)).T)).T[0]
    return r*M * (1 - np.dot(A, M) / k)

def getRandomMatrices(N):
    # Definition of parameters
    
    # M0 => intial population counts
    # initial abundances of the species were drawn from a uniform distribution ranging between 10 and 100.
    x0 = []
    for i in range(0,N):
        x0.append(random.randint(10,100))
    x0 = np.array(x0)
    
    # B => growth rates
    # Growth rates were assigned to each species from a uniform distribution from 0 to 1 so that all species were capable of positive growth
    r = []
    for i in range(0,N):
        r.append(random.uniform(0.01, 1.0))
    r = np.array(r)
    
    # CARRYING CAPACITY beta distribution, but also normal distribution.
    # scaled to range from 1 to 100
    k = []
    for i in range(0,N):
        k.append(random.randint(1,100))
    k = np.array(k)

    # A => interaction matrix, i=j should only be negative
    A = []
    return A, r, k, x0

def getK(N):
    k = []
    for i in range(0,N):
        k.append(random.randint(1,100))
    k = np.array(k)
    return k

def getR(N):
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

def getEnvCommunities(A):
    G = nx.from_pandas_adjacency(A)

def step1_createLoadA(workdir, base_graph_path):
    sim_A = pd.read_csv(base_graph_path, delimiter=",", header=0, index_col=0)
    sim_A.index = sim_A.columns
    
    #print ("generating new values FOR simulation matrix")
    N = sim_A.shape[0]
    multM = np.random.uniform(low=.5, high=1.0, size=(N,N))#np.random.rand(N,N) # low=.5, high=1.0, size=(N,N)
    
    # get random matrix of [-1,1]. Symmetrical means no predation
    randInts = np.asmatrix(np.reshape(np.random.choice([-1,1], N*N), (N,N) )) #[-1,1]
    randInts_symm = np.triu(randInts, k=1) + np.triu(randInts, k=1).T
    
    new_sim_A = sim_A * multM * randInts #randInts_symm
    new_sim_A = new_sim_A.replace([-0], 0)
    new_sim_A = round(new_sim_A,1)
    
    new_sim_A.index = pd.Index([k for k in range(1,N+1)]) 
    new_sim_A.columns = pd.Index([k for k in range(1,N+1)]) 
    new_sim_A = fillDiagonal(new_sim_A, value = 1.0)

    new_sim_A.index = new_sim_A.columns
    return new_sim_A, sim_A
    
def step2_createLoadTemps(thisA, existing_tpos, existing_tnegs, choicetype): #"cluster", "random"
    posindeces = []
    negindeces = []
    
    if existing_tpos and existing_tnegs:
        
        tpos = list(pd.read_csv(existing_tpos, delimiter=",", header=0, index_col=0).index)
        tneg = list(pd.read_csv(existing_tnegs, delimiter=",", header=0, index_col=0).index)
        
        for pos in tpos:
            if pos in list(thisA.index):
                posindeces.append(list(thisA.index).index(pos))
        for neg in tneg:
            if neg in list(thisA.index):
                negindeces.append(list(thisA.index).index(neg))
    else:
        
        if choicetype=="cluster":
            G = nx.from_pandas_adjacency(thisA)
            mod, comms = modularity(G)
            
            coldwarm = "warm"
            for comm in comms:
                if len(comm) == 2:
                    n_choices = 2
                else:
                    n_choices = int(len(comm)/1.2)
                indeces = np.random.choice(list(comm), n_choices, replace=False) # select (half) of the size of the community to be associated to temperature
                if coldwarm == "warm":
                    for ind in indeces:
                        posindeces.append(ind)
                    coldwarm = "cold"
                else:
                    for ind in indeces:
                        negindeces.append(ind)
                    coldwarm = "warm"
        if choicetype=="random":
            n_choices = int(len(list(thisA.columns))*(2/3))
            indeces = np.random.choice(list(thisA.columns), n_choices, replace=False)
            posindeces = indeces[:int(n_choices/2)]
            negindeces = indeces[int(n_choices/2):]

    return posindeces, negindeces

def filtMatrices(abunds, new_sim_A, base_sim_A):
    
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

workdir = snakemake.config['workdir']
nsimulations = snakemake.config['nsimulations']
#outputdir = snakemake.output[0]
base_graph_dir = snakemake.input[0]

#workdir = "C:/Users/vdinkel/Desktop/Data/4_SyntheticGenerators/data/"

# example simulation
n_samples = 100
t = np.linspace(0.0, 100.0,  100)

nrun = 0
while nsimulations >0:
    print ("--- nsimulations: ",nsimulations)
    new_sim_A, base_sim_A = step1_createLoadA(workdir, base_graph_dir)
    N = new_sim_A.shape[0]
    A, r, k, x0 = getRandomMatrices(N)
    A = new_sim_A
    
    k = getK(N)
    r = getR(N)
    X, infodict = integrate.odeint(dX_dt_gLV, x0, t, full_output=True)
    #plt.ylim(0, 1000)
    #plt.plot(X)
    #plt.show()
    
    testV = np.matrix(abs(X)).flatten()
    testV = np.array(testV)[0]
    
    if max(testV) > 10000 or sum(X[99]) == 0:
        print ("-- FAILED SIM")
    else:
        print ("-- GOOD SIM")
        plt.plot(X)
        plt.title("Single Simulation")
        plt.xlabel("Time")
        plt.ylabel("Abundance")
        plt.savefig(snakemake.output[nrun].replace(".csv","_sim_.jpg"))
        plt.clf()

        abunds = []
        for j in range(0, n_samples):
            k = getK(N)
            r = getR(N)
            X, infodict = integrate.odeint(dX_dt_gLV, x0, t, full_output=1)
            abunds.append(X[99])
        
        #import pdb; pdb.set_trace()
        if (max(np.array(abs(np.matrix(abunds).flatten()))[0]) > 10000) or sum(pd.DataFrame(abunds).T.sum()>1) <20:
            print ("-- FAILED ABUNDS") 
        else:
            print ("-- GOOD ABUNDS") 
            plt.plot(abunds)
            plt.title("Simulated Attractor Abundances")
            plt.xlabel("Simulation")
            plt.ylabel("Attracted Abundance")
            plt.savefig(snakemake.output[nrun].replace(".csv","_abunds_.jpg"))
            plt.clf()
    
            abunds = pd.DataFrame(abunds)
            abunds = abunds.where(abunds >= 1, 0)
    
            # STEP 4: nullsummen entfernen
            filt_abunds, filt_base_A, filt_new_sim_A, badcols = filtMatrices(abunds, fillDiagonal(new_sim_A,0) , base_sim_A) #fillDiagonal(new_sim_A, 0)
            print ("-- REMOVED ", str(len(badcols)), " Indeces")
    
            # STEP 5: alles speichern
            pd.DataFrame.to_csv(abunds, snakemake.output[nrun], index=True, sep=",")
            
            pd.DataFrame.to_csv(pd.DataFrame(new_sim_A), snakemake.output[nrun].replace(".csv","_new_sim_A.csv"))
            pd.DataFrame.to_csv(pd.DataFrame(base_sim_A), snakemake.output[nrun].replace(".csv","_base_sim_A.csv"))
            
            pd.DataFrame.to_csv(pd.DataFrame(filt_abunds), snakemake.output[nrun].replace(".csv","_filt_abunds.csv"))
            pd.DataFrame.to_csv(pd.DataFrame(filt_base_A), snakemake.output[nrun].replace(".csv","_filt_base_A.csv"))
            pd.DataFrame.to_csv(pd.DataFrame(filt_new_sim_A), snakemake.output[nrun].replace(".csv","_filt_new_sim_A.csv"))
            
            #ESABO
            matrix = binNorm(filt_abunds)
            esaboM = corrESABOM(matrix, 1.3)
            absEsabo = abs(esaboM.round().astype(int))
            absEsabo = absEsabo.where(absEsabo == 0, 1)
            #import pdb; pdb.set_trace()
            
            #simV = fillDiagonal(filt_base_A, 0)
            #simV = np.matrix(simV).flatten()
            #simV = np.array(simV)[0]
            #esaboV = np.array(np.matrix(esaboM).flatten())[0]
            #corr, pVal = stats.pearsonr(simV, esaboV)
            
            pd.DataFrame.to_csv(absEsabo, snakemake.output[nrun].replace("/abundances/","/networks/").replace(".csv","_esabo.csv"))
            
            nsimulations -= 1
            nrun += 1