# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:29:44 2023

@author: vdinkel
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

def msum(A):
    return int(np.sum(np.sum(A)))

def AND(A, B):
    return B.where(A == 1, 0)

def fillDiagonal(dataFrame, value = 0):
    i = 0
    for cols in dataFrame.columns:
        dataFrame.at[dataFrame.index[i], cols] = value
        i += 1
    return dataFrame

def getOwnTPs(trueA, predA, threshold = 0.0):
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
    if tp+fp+tn+fn != 0:
        acc = (tp+tn) / (tp+fp+tn+fn)
    else:
        acc = 0
    return {"tp": tp, "fp":fp, "tn": tn, "fn": fn, "tpr": tpr, "fpr": fpr, "acc": acc}

def getPosNegAs(thisA):
    A_neg = thisA.where( thisA < 0, 0)
    A_neg = A_neg.where( A_neg == 0, 1)
    A_pos = thisA.where( thisA > 0, 0)
    A_pos = A_pos.where( A_pos == 0, 1)
    return A_pos, A_neg

def getMultiplexThresh(multA, thresh):
    multA = multA.where(multA >= thresh, 0)
    multA = multA.where(multA == 0, 1)
    return multA

def getMultiplexReverseThresh(multA, thresh):
    #import pdb; pdb.set_trace()
    multA = multA.where(multA < thresh, 0)
    multA = multA.where(multA == 0, 1)
    return multA

def processOutputs():
    # this function aggregates the simulation runs into single quality statistics e.g. accuracy, true positive rates etc.
    # in the pipeline output, only positive predictive value (PPV) and the rate of true positives to false positives (TP/FP) are relevant
    # this aggregation also includes analysis of positive and negative links and reverse runs.
    # for the sake of completeness, these additional aggregations are kept in with increased computation time as tradeoff
    outputs = []
    allRuns = {}
    allAccs = {}
    allTprs = {}
    allFprs = {}
    allPpvs = {}
    allCsis = {}
    allReverseRuns = {}
    allReverseAccs = {}
    allReverseTprs = {}
    allReverseFprs = {}
    allReversePpvs = {}
    allReverseCsis = {}
    for meth in methods_list:
        allRuns[meth] = []
        allAccs[meth] = []
        allTprs[meth] = []
        allFprs[meth] = []
        allPpvs[meth] = []
        allCsis[meth] = []
        # REVERSE
        allReverseRuns[meth] = []
        allReverseAccs[meth] = []
        allReverseTprs[meth] = []
        allReverseFprs[meth] = []
        allReversePpvs[meth] = []
        allReverseCsis[meth] = []
    
    for i in range(0, nruns):
        print (i)
        inputFileDir = basedir+"glv_"+nettype+"_"+str(i)+"_filt_base_A.csv"
    
        # LOAD BASE SIMULATION NETWORKS
        sim_A = pd.read_csv(inputFileDir, delimiter=",", header=0, index_col=0)
        A_vals = pd.read_csv(inputFileDir.replace("filt_base_A","filt_new_sim_A"), delimiter=",", header=0, index_col=0)
        
        # LOAD INFERRED NETWORKS
        inputNetsDir = inputFileDir.replace("abundances","networks").replace("filt_base_A.csv","")
        propr_A = pd.read_csv(inputNetsDir+"propr.csv", delimiter=",", header=0, index_col=0)
        ccrepe_A = pd.read_csv(inputNetsDir+"ccrepe.csv", delimiter=",", header=0, index_col=0)
        sparcc_A = pd.read_csv(inputNetsDir+"sparcc.csv", delimiter=",", header=0, index_col=0)
        spearman_A = pd.read_csv(inputNetsDir+"spearman.csv", delimiter=",", header=0, index_col=0)
        spieceasi_A = pd.read_csv(inputNetsDir+"spieceasi.csv", delimiter=",", header=0, index_col=0)
        esabo_A =  pd.read_csv(inputNetsDir+"esabo.csv", delimiter=",", header=0, index_col=0)
        ecocopula_A =  pd.read_csv(inputNetsDir+"ecocopula.csv", delimiter=",", header=0, index_col=0)
        
        # set the diagonal to 0 as these are no interactions
        A_vals = fillDiagonal(A_vals, 0)
        A_vals_pos, A_vals_neg = getPosNegAs(A_vals)
        A_vals_pos.index = propr_A.index
        A_vals_pos.columns = propr_A.columns
        A_vals_neg.index = propr_A.index
        A_vals_neg.columns = propr_A.columns
        
        # clean up row and column names
        sim_A.index = propr_A.index
        sim_A.columns = propr_A.columns
        
        esabo_A.index = propr_A.index
        esabo_A.columns = propr_A.columns
        esabo_A = esabo_A.where(abs(esabo_A) == 0, 1)
        
        abs_sim_A = sim_A.where(abs(sim_A) == 0, 1)
 
        # generate consensus network
        n_multiplex = 7
        multiplex_A = propr_A + spieceasi_A + ccrepe_A + sparcc_A + spearman_A + esabo_A + ecocopula_A
        multiplex_A = multiplex_A / n_multiplex
        w_multiplex_A = propr_A + spieceasi_A + ccrepe_A + sparcc_A + spearman_A + esabo_A + ecocopula_A
        w_multiplex_A = w_multiplex_A/np.max(w_multiplex_A.values)
        
        ownTPs = []
        w_ownTPs = []
        # go through the consensus threshold 0.0 to 1.0 with 0.1 steps and asses TPs and FPs
        for k in range(0, 11):
            ownTPs.append(getOwnTPs(abs_sim_A, getMultiplexThresh(multiplex_A, k/10)))
        for k in range(0, 11):
            w_ownTPs.append(getOwnTPs(abs_sim_A, getMultiplexThresh(w_multiplex_A, k/10)))
    
        # get consensus network at configured (optimal) threshold of 0.5
        newBest = getMultiplexThresh(w_multiplex_A, 0.5)
        
        allmethods = {"propr_A": propr_A, "spieceasi_A": spieceasi_A, "ccrepe_A": ccrepe_A, "sparcc_A": sparcc_A, "spearman_A": spearman_A, "esabo_A": esabo_A, "ecocopula_A": ecocopula_A, "mask_A": getMultiplexThresh(multiplex_A, 0.0)} #

        fnrs = {}
        tpfn_ratios = {}
        tpfp_ratios = {}
        

        fnrs_reverse = {}
        tpfn_ratios_reverse = {}
        tpfp_ratios_reverse = {}
        
        # go through all methods and consensus threshold 0.0 to 1.0 and get inference quality metrics
        for key in allmethods.keys():
            method = allmethods[key]
            
            thisTPsL = []
            thisFPsL = []
            thisFNsL = []
            thisACCsL = []
            thisTPRsL = []
            thisFPRsL = []
            thisPPVsL = []
            thisCSIsL = []
            thisReverseTPsL = []
            thisReverseFPsL = []
            thisReverseFNsL = []
            thisReverseACCsL = []
            thisReverseTPRsL = []
            thisReverseFPRsL = []
            thisReversePPVsL = []
            thisReverseCSIsL = []
            for j in range(0,11):
                mask_net_t = getMultiplexThresh(multiplex_A, j/10)
                mask_reverse_net_t = getMultiplexReverseThresh(multiplex_A, j/10)
                thisNet = AND(method,mask_net_t)
                thisReverseNet = AND(method,mask_reverse_net_t)
                thisTPs = getOwnTPs(abs_sim_A, thisNet)
                thisReverseTPs = getOwnTPs(abs_sim_A, thisReverseNet)
                thisTPsL.append(thisTPs['tp'])
                thisFPsL.append(thisTPs['fp'])
                thisFNsL.append(thisTPs['fn'])
                thisACCsL.append(thisTPs['acc'])
                #thisTPRsL.append(thisTPs['tpr'])
                try:
                    thisTPRsL.append((thisTPs['tp']/thisTPs['fp'])/(thisTPs['fn']/thisTPs['tn']))
                except:
                    thisTPRsL.append(0)
                thisFPRsL.append(thisTPs['fpr'])
                if (thisTPs['tp']+thisTPs['fp']) > 0:
                    thisPPVsL.append(thisTPs['tp'] / (thisTPs['tp']+thisTPs['fp']))
                else:
                    thisPPVsL.append(0)
                if (thisTPs['tp']+thisTPs['fn']+thisTPs['fp']) > 0:
                    thisCSIsL.append(thisTPs['tp'] / (thisTPs['tp']+thisTPs['fn']+thisTPs['fp']))
                else:
                    thisCSIsL.append(0)
                
                thisReverseTPsL.append(thisReverseTPs['tp'])
                thisReverseFPsL.append(thisReverseTPs['fp'])
                thisReverseFNsL.append(thisReverseTPs['fn'])
                thisReverseACCsL.append(thisReverseTPs['acc'])
                try:
                    thisReverseTPRsL.append((thisReverseTPs['tp']/thisReverseTPs['fp'])/(thisReverseTPs['fn']/thisReverseTPs['tn']))
                except:
                    thisReverseTPRsL.append(0)
                thisReverseFPRsL.append(thisReverseTPs['fpr'])
                if (thisReverseTPs['tp']+thisReverseTPs['fp']) > 0:
                    thisReversePPVsL.append(thisReverseTPs['tp'] / (thisReverseTPs['tp']+thisReverseTPs['fp']))
                else:
                    thisReversePPVsL.append(0)
                if (thisReverseTPs['tp']+thisReverseTPs['fn']+thisReverseTPs['fp']) > 0:
                    thisReverseCSIsL.append(thisReverseTPs['tp'] / (thisReverseTPs['tp']+thisReverseTPs['fn']+thisReverseTPs['fp']))
                else:
                    thisReverseCSIsL.append(0)
            
            fnrs[key] = [1.0-k for k in thisTPRsL]
            tpfn_ratios[key] = np.array(thisTPsL) / np.array(thisFNsL)
            
            thisFPsL = [1 if x == 0 else x for x in thisFPsL]
            tpfp_ratio = np.array(thisTPsL)/np.array(thisFPsL)
            tpfp_ratios[key] = tpfp_ratio
            
            allRuns[key.replace("_A","")].append(tpfp_ratio)
            allAccs[key.replace("_A","")].append(thisACCsL)
            allTprs[key.replace("_A","")].append(thisTPRsL)
            allFprs[key.replace("_A","")].append(thisFPRsL)
            allPpvs[key.replace("_A","")].append(thisPPVsL)
            allCsis[key.replace("_A","")].append(thisCSIsL)
            
            fnrs_reverse[key] = [1.0-k for k in thisReverseTPRsL]
            tpfn_ratios_reverse[key] = np.array(thisReverseTPsL) / np.array(thisReverseFNsL)
            thisReverseFPsL = [1 if x == 0 else x for x in thisReverseFPsL]
            tpfp_ratio_reverse = np.array(thisReverseTPsL)/np.array(thisReverseFPsL)
            tpfp_ratios_reverse[key] = tpfp_ratio_reverse
            
            allReverseRuns[key.replace("_A","")].append(tpfp_ratio_reverse)
            allReverseAccs[key.replace("_A","")].append(thisReverseACCsL)
            allReverseTprs[key.replace("_A","")].append(thisReverseTPRsL)
            allReverseFprs[key.replace("_A","")].append(thisReverseFPRsL)
            allReversePpvs[key.replace("_A","")].append(thisReversePPVsL)
            allReverseCsis[key.replace("_A","")].append(thisReverseCSIsL)
        
        A_vals.index = allmethods['spieceasi_A'].index
        A_vals.columns = allmethods['spieceasi_A'].columns
        G = nx.from_pandas_adjacency(A_vals, create_using=nx.DiGraph)
        
        inttypes = {}
        for key in allmethods.keys():
            inttypes[key] = {"pos": 0, "posneg": 0, "neg": 0, "edges": []}
        inttypes["base"] = {"pos": 0, "posneg": 0, "neg": 0}
        for edge in G.edges:
            inttype = ""
            if A_vals[edge[1]][edge[0]] > 0 and A_vals[edge[0]][edge[1]] > 0:
                # both positive
                inttype = "pos"  
            if (A_vals[edge[1]][edge[0]] > 0 and A_vals[edge[0]][edge[1]] < 0) or (A_vals[edge[1]][edge[0]] < 0 and A_vals[edge[0]][edge[1]] > 0):
                # pos/neg interaction
                inttype = "posneg"
            if A_vals[edge[1]][edge[0]] < 0 and A_vals[edge[0]][edge[1]] < 0:
                # both negative
                inttype = "neg"
            inttypes["base"][inttype] += 1
            for key in allmethods.keys():
                if allmethods[key][edge[0]][edge[1]] > 0:
                    inttypes[key][inttype] += 1
                    inttypes[key]["edges"].append([A_vals[edge[1]][edge[0]], A_vals[edge[1]][edge[0]]]) 
        
        outputs.append([abs_sim_A, allmethods])
    
    # save the aggregated inference metrics into pickle files for future loading
    with open(benchmarkdir+'allRuns_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allRuns, f)
    with open(benchmarkdir+'allAccs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allAccs, f)
    with open(benchmarkdir+'allTprs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allTprs, f)
    with open(benchmarkdir+'allFprs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allFprs, f)
    with open(benchmarkdir+'allPpvs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allPpvs, f)
    with open(benchmarkdir+'allCsis_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allCsis, f)
    with open(benchmarkdir+'allReverseRuns_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allReverseRuns, f)
    with open(benchmarkdir+'allReverseAccs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allReverseAccs, f)
    with open(benchmarkdir+'allReverseTprs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allReverseTprs, f)
    with open(benchmarkdir+'allReverseFprs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allReverseFprs, f)
    with open(benchmarkdir+'allReversePpvs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allReversePpvs, f)
    with open(benchmarkdir+'allReverseCsis_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(allReverseCsis, f)
    with open(benchmarkdir+'allOutputs_'+str(nruns-1)+'.pkl', 'wb') as f: 
        pickle.dump(outputs, f)
        
    return outputs

def loadProcessedOutputs(benchmarkdir):
    # loads the processed pickle files (aggregated inference quality metrics)
    file_allRuns = open(benchmarkdir+'allRuns_'+str(nruns-1)+'.pkl', 'rb')
    allRuns = pickle.load(file_allRuns)
    file_allRuns.close()
    
    file_allAccs = open(benchmarkdir+'allAccs_'+str(nruns-1)+'.pkl', 'rb')
    allAccs = pickle.load(file_allAccs)
    file_allAccs.close()
    
    file_allTprs = open(benchmarkdir+'allTprs_'+str(nruns-1)+'.pkl', 'rb')
    allTprs = pickle.load(file_allTprs)
    file_allTprs.close()
    
    file_allFprs = open(benchmarkdir+'allFprs_'+str(nruns-1)+'.pkl', 'rb')
    allFprs = pickle.load(file_allFprs)
    file_allFprs.close()
    
    file_allPpvs = open(benchmarkdir+'allPpvs_'+str(nruns-1)+'.pkl', 'rb')
    allPpvs = pickle.load(file_allPpvs)
    file_allPpvs.close()
    
    file_allCsis = open(benchmarkdir+'allCsis_'+str(nruns-1)+'.pkl', 'rb')
    allCsis = pickle.load(file_allCsis)
    file_allCsis.close()
    
    file_allReverseRuns = open(benchmarkdir+'allReverseRuns_'+str(nruns-1)+'.pkl', 'rb')
    allReverseRuns = pickle.load(file_allReverseRuns)
    file_allReverseRuns.close()
    
    file_allReverseAccs = open(benchmarkdir+'allReverseAccs_'+str(nruns-1)+'.pkl', 'rb')
    allReverseAccs = pickle.load(file_allReverseAccs)
    file_allReverseAccs.close()
    
    file_allReverseTprs = open(benchmarkdir+'allReverseTprs_'+str(nruns-1)+'.pkl', 'rb')
    allReverseTprs = pickle.load(file_allReverseTprs)
    file_allReverseTprs.close()
    
    file_allReverseFprs = open(benchmarkdir+'allReverseFprs_'+str(nruns-1)+'.pkl', 'rb')
    allReverseFprs = pickle.load(file_allReverseFprs)
    file_allReverseFprs.close()
    
    file_allReversePpvs = open(benchmarkdir+'allReversePpvs_'+str(nruns-1)+'.pkl', 'rb')
    allReversePpvs = pickle.load(file_allReversePpvs)
    file_allReversePpvs.close()
    
    file_allReverseCsis = open(benchmarkdir+'allReverseCsis_'+str(nruns-1)+'.pkl', 'rb')
    allReverseCsis = pickle.load(file_allReverseCsis)
    file_allReverseCsis.close()
    
    return allRuns, allAccs, allTprs, allFprs, allPpvs, allCsis, allReverseRuns, allReverseAccs, allReverseTprs, allReverseFprs, allReversePpvs, allReverseCsis
      

def plotBenchmark(axs, row, runIndeces, center, runDict, method, y_label, color, isReverse):
    
    base = []
    maxdatas = []
    maxindeces = []
    
    # processed data contains lists of runs and increasing thresholds: x[0] = run1[0], run1[1], run1[2], ...
    # Because we want to plot results of one threshold of all runs, we have to unfold so that x[0] = run1[0], run2[0], etc...
    
    # alreadyGoodIndeces are indeces of runs where TP > FP at initial inference
    alreadyGoodIndeces = list(set([k for k in range(0,nruns)]) - set(runIndeces[method]['badRunIndeces']))
    # unswitchedIndeces are indeces of runs where FP > TP at all consensus thresholds
    unswitchedIndeces = list(set(runIndeces[method]['badRunIndeces']) - set(runIndeces[method]['switchRunIndeces']))
    # switchedIndeces are indeces of runs where initial TP < FP and with consensus network became TP>FP
    switchedIndeces = runIndeces[method]['switchRunIndeces']

    indecesList = [alreadyGoodIndeces, switchedIndeces, unswitchedIndeces]
    origDict = runDict[method]
    titles = ["Initial Inference with TP>FP ("+str(int(100*(len(alreadyGoodIndeces)/nruns)))+"%)", "Switched inital TP<FP to TP>FP ("+str(int(100*(len(switchedIndeces)/1000)))+"%)", "Bad Initial Inference (TP<FP): "+str(len(unswitchedIndeces))]
    
    # this is the separation of results into alreadyGoodIndeces and switchedIndeces (unswitchedIndeces are disregarded)
    j = 0
    for onlyIndeces in indecesList[:-1]:
        
        data = []
        print (method, " ", len(onlyIndeces))
        runDict[method] = origDict
        runDict[method] = [runDict[method][x] for x in onlyIndeces]
        
        if center:
            i = 0
            for run in runDict[method]:

                centRun = list(np.array(run) - run[0])
                centRun[0] = run[0]
                runDict[method][i] = centRun
                i += 1
        base.append([k[0] for k in runDict[method]])    
        

        for i in range(1, 11):
            # UNFOLD DATA
            data.append([k[i] for k in runDict[method]])

        maxdatas.append(np.array([k[0] for k in runDict[method]]) + [max(k[1:]) for k in runDict[method]])
        maxindeces.append([list(k[1:]).index(max(k[1:])) for k in runDict[method]])
        meanSeries = np.array([np.mean(k) for k in data])
        stdSeries = np.array([np.std(k) for k in data])
        
        stdMaxSeries = []
        stdMinSeries = []
        for dat in data:
            posonly = [k for k in dat if k>=0]
            negonly = [k for k in dat if k<=0]
            stdMaxSeries.append(np.std(posonly))
            stdMinSeries.append(np.std(negonly))
        
        stdMaxSeries = [k if not np.isnan(k) else 0 for k in stdMaxSeries]
        stdMinSeries = [k if not np.isnan(k) else 0 for k in stdMinSeries]

        axs[row, j].plot(meanSeries, linewidth = 3, color="black", alpha=0.5, label="mean")
        
        y1 = meanSeries + abs(np.array(stdMaxSeries))
        y2 = meanSeries - abs(np.array(stdMinSeries))
 
        axs[row, j].fill_between([i for i in range(0, 10)], y1, meanSeries, where=y1 >= meanSeries, facecolor=(0.182, 0.750, 0.172, .2), interpolate=True, label="std of positives")
        axs[row, j].fill_between([i for i in range(0, 10)], y2, meanSeries, where=y2 <= meanSeries, facecolor=color, interpolate=True, label="std of negatives")

        axs[row, j].axvline(x = 0, color = 'gray', linestyle = "solid", linewidth=0.6)
        axs[row, j].axhline(y = 0, color = 'gray', linestyle = 'solid', linewidth=0.6)
        
        if row == 0:
            axs[row, j].set_title(titles[j])
            axs[row, j].set_xticklabels(["0.0", "Initial\n Inference", "0.2", "0.4", "0.6", "0.8"])
        if j == 0:
            axs[row, j].set(ylabel=y_label)
        if row == 1:
            axs[row, j].set(xlabel="Ensemble Network Threshold")
            
        j += 1
    axs[0, 1].legend()   
    
    init_tpfp = (int(np.mean([item for sublist in base for item in sublist])*100))/100
    init_tpfp_good = (int((len(alreadyGoodIndeces)/nruns)*100))/100
    masked_tpfp = (int((np.mean([item for sublist in maxdatas for item in sublist]) + np.mean([item for sublist in base for item in sublist]))*100))/100
    masked_tpfp_good = (int(((len(alreadyGoodIndeces) + len(switchedIndeces))/nruns)*100))/100
    opt_masked_indeces = (int((np.mean([item for sublist in maxindeces for item in sublist]))*100))/100
    print ("Init. TP/FP: ", init_tpfp)
    print ("% Init. TP/FP>1: ", init_tpfp_good)
    print ("Masked. TP/FP: ", masked_tpfp)
    print ("% Masked. TP/FP>1: ", masked_tpfp_good)
    print ("Opt. Masked Indeces: ", opt_masked_indeces)
    return {"initial tp/fp": init_tpfp, "proportion of initial tp>fp": init_tpfp_good, "cn tp/fp": masked_tpfp, "proportion of cn tp>fp":masked_tpfp_good, "optimal en threshold":opt_masked_indeces}

# load snakemake configuration
nruns = snakemake.config["nsimulations"]
nettype = snakemake.config["nettype"]
seed = snakemake.config["seed"]
workdir = snakemake.config["workdir"]+"outputs/"+str(seed)+"/"

basedir = workdir+"abundances/"
netsdir = workdir+"networks/" 
benchmarkdir = workdir+"benchmark/" 
plotdir = benchmarkdir
supplementsdir = benchmarkdir
methods_list = ["spieceasi", "esabo", "ccrepe", "sparcc", "propr", "spearman", "ecocopula", "mask"] #consensus network was once called mask

inputFileDir = basedir+"glv_"+nettype+"_"+str(nruns-1)+"_filt_base_A.csv"

if True: # set here to False if the outputs were already processed once
    outputs = processOutputs()

print ("Loading processed files")
allRuns, allAccs, allTprs, allFprs, allPpvs, allCsis, allReverseRuns, allReverseAccs, allReverseTprs, allReverseFprs, allReversePpvs, allReverseCsis = loadProcessedOutputs(benchmarkdir)  
print("... Done")

sns.set(rc={'figure.figsize':(11.7,8.27)})


# 1. Find indeces where switch TP>FP happens
# the result of this function is not part of the pipeline but kept for completeness
p = 0
runIndeces = {}

for meth in allRuns.keys():

    plt.figure(p)
    badRunIndeces = [index for (index, item) in enumerate([allRuns[meth][i][0]<1 for i in range(0,nruns)]) if item]
    switchRunIndeces = [index for (index, item) in enumerate([allRuns[meth][i][0]<1 and True in (allRuns[meth][i] > 1) for i in range(0,nruns)]) if item]
    switchIndeces = [(list(allRuns[meth][i] < 1)).index(False) for i in switchRunIndeces] # find indeces where ratio is smaller 1
    runIndeces[meth] = {"badRunIndeces": badRunIndeces, "switchRunIndeces": switchRunIndeces}
    
    counts, bins = np.histogram(switchIndeces)
    plt.stairs(counts, bins, fill=True)
    totalBadRunPercentage = int(np.round(len(badRunIndeces) / nruns, 2) * 100)
    try:
        switchedPercentage = int(np.round(len(switchRunIndeces) / len(badRunIndeces), 2)*100)
    except:
        switchedPercentage = 0
    improvedBadRunPercentage = int(np.round((len(badRunIndeces) - len(switchRunIndeces)) / nruns, 2)*100)
    plt.title(meth.upper()+"\nTP<FP percentage: "+str(totalBadRunPercentage)+"% ~> "+str(improvedBadRunPercentage)+" %\nTP>FP switch percentage: "+str(switchedPercentage)+" %")
    p += 1

# plot the TP/FPs and PPVs 
method_statistics = {}
method_ppv_statistics = {}
for method in methods_list:
    
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    fig, axs = plt.subplots(2, 2, sharex = True, sharey = "row", constrained_layout=True)
    title = method
    print (method)
    if method == "mask":
        title = "Ensemble Network"
    if method == "spieceasi":
        title = "spiec-easi"
    fig.suptitle(title.upper()+"\nRelative Precision Improvement", fontsize=25)
    
    doCenter = True # centering means using increase/decrease of quality improvement centered (in relation) to its initial inference
    tpfps = plotBenchmark(axs,0, runIndeces, doCenter, allRuns, method, 'TP/FP', (.4, .6, .8, .5), False)
    ppvs = plotBenchmark(axs, 1, runIndeces, doCenter, allPpvs, method, 'PPV', (.4, .6, .8, .5), False)

    p += 1

    method_statistics[method] = tpfps
    method_ppv_statistics[method] = ppvs
    
    if method == "mask":
        method = "cn"
    fig.savefig(plotdir+'/benchmark_'+method+'_initial_cn.png', format='png', dpi=1200)


pd_statistics = pd.DataFrame(method_statistics).T
pd_statistics.to_csv(supplementsdir+"/method_synthetic_benchmark.csv", sep=',', index=True, encoding='utf-8')