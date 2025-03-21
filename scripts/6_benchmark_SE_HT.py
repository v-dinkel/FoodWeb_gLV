# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:07:41 2023

@author: vdinkel
"""
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

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

def msum(A):
    return int(np.sum(np.sum(A)))

def AND(A, B):
    return B.where(A == 1, 0)

def OR(A, B):
    C = abs(A) + abs(B)
    return C.replace(2, 1)

def fillDiagonal(dataFrame, value = 0):
    i = 0
    for cols in dataFrame.columns:
        dataFrame.at[dataFrame.index[i], cols] = value
        i += 1
    return dataFrame

def getMultiplexThresh(multA, thresh):
    multA = multA.where(multA >= thresh, 0)
    multA = multA.where(multA == 0, 1)
    return multA

# set the paths according to the seed which you want to compare
simpath = "/home/viktor/project_migration/FoodWeb_gLV/outputs/100/abundances/"
path = "/home/viktor/project_migration/FoodWeb_gLV/outputs/100/networks/"
benchpath = "/home/viktor/project_migration/FoodWeb_gLV/outputs/100/benchmark/"

# don't change these filebases
filebase = "glv_cluster_<i>_<method>.csv"
simbase = "glv_cluster_<i>_filt_base_A.csv"

tps_diffs = []
fps_diffs = []
for i in range(0,10):
    
    # load the simulation adjacency matrix (ground truth)
    file = filebase.replace("<i>",str(i))
    simfile = simbase.replace("<i>",str(i))
    sim_A = pd.read_csv(simpath+simfile, delimiter=",", header=0, index_col=0)
    
    # load the networks inferred by different inference methods
    propr = pd.read_csv(path+file.replace("<method>", "propr"), delimiter=",", header=0, index_col=0)
    ccrepe = pd.read_csv(path+file.replace("<method>", "ccrepe"), delimiter=",", header=0, index_col=0)
    spieceasi = pd.read_csv(path+file.replace("<method>", "spieceasi"), delimiter=",", header=0, index_col=0)
    sparcc = pd.read_csv(path+file.replace("<method>", "sparcc"), delimiter=",", header=0, index_col=0)
    esabo = pd.read_csv(path+file.replace("<method>", "esabo"), delimiter=",", header=0, index_col=0)
    ecocopula = pd.read_csv(path+file.replace("<method>", "ecocopula"), delimiter=",", header=0, index_col=0)
    spearman = pd.read_csv(path+file.replace("<method>", "spearman"), delimiter=",", header=0, index_col=0)
    
    methods_df = {"propr": propr, "ccrepe": ccrepe, "spieceasi": spieceasi,"esabo": esabo, "ecocopula": ecocopula, "sparcc": sparcc, "spearman": spearman}
    
    # create the ensemble network, called multiplex_A here
    multiplex_A =  methods_df['propr'] + methods_df['spieceasi'] + methods_df['esabo'] + methods_df['sparcc'] + methods_df['spearman'] + methods_df['ecocopula'] + methods_df['ccrepe']
    nmult = 7
    multiplex_A = multiplex_A / nmult
    
    # define base network and create intersection with ensemble to return consensus network
    method = "spieceasi"
    baseNet = methods_df[method]
    multNet = AND(methods_df[method], getMultiplexThresh(multiplex_A, 0.5))
    
    # read the corresponding weighted spiec easi network
    spieceasi_HT = pd.read_csv(path+file.replace("<method>", "spieceasi_weighted"), delimiter=",", header=0, index_col=0)
    
    # find threshold for association strength in spiec easi network, which comes closest to the edge count of the consensus network
    for j in range(0, 100) :
        t = j/100
        ht_net = getMultiplexThresh(spieceasi_HT, t)
        if msum(ht_net) <= msum(multNet):
            break
    
    # compare consensus network and spiec easi network with ground truth network to get true positives and false positives
    ret_mult = getOwnTPs(sim_A, multNet, 0.0)
    ret_ht = getOwnTPs(sim_A, ht_net, 0.0)
    
    # store difference in TPs and FPs
    tps_diffs.append(ret_mult['tp'] - ret_ht['tp'])
    fps_diffs.append(ret_mult['fp'] - ret_ht['fp'])

# plot
df_tp = pd.DataFrame(np.array([tps_diffs, fps_diffs])).T
df_tp.boxplot()
plt.title("Consensus Network compared to SpiecEasi")
plt.xticks([1,2], ["TPs", "FPs"])
plt.ylabel("Difference as CN - SE")

# save plot
plt.savefig(benchpath+"compare_CN_SE_HT.png")