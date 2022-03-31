# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:14:11 2021

@author: leokt
"""
import pandas as pd
import numpy as np
import math
import requests as rq
from io import BytesIO


# =============================================================================
# module for generating dot-score rankings
# =============================================================================
def dotscores(kin_vals, clus_vals, kin_names, clus_names):
    y_scores = list(range(len(clus_vals)))
    for u in range(len(clus_vals)):
        x_scores = list(range(len(kin_vals)))
        for p in range(len(kin_vals)):
            rank = np.dot(kin_vals[p], clus_vals[u])
            x_scores[p] = rank
        y_scores[u] = x_scores
    xydot_df = pd.DataFrame(np.transpose(y_scores), columns = clus_names, index = kin_names)
    return xydot_df

# =============================================================================
# interpolating scores from KinMap nearest neighbor match
# kinmatrix refers to the sugiyama-derived matrix of interest (just one)
# =============================================================================
def interpolate(kinmatrix):
    url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/interpolator_kinases.xlsx?raw=true"
    data = rq.get(url).content
    interpoltools = pd.read_excel(BytesIO(data), None)
    kinmatrix.reset_index(inplace=True)
    kinmatrix = kinmatrix.rename(columns = {'index':'Human Kinase'})
    merge1 = kinmatrix.merge(interpoltools['match 355'], on='Human Kinase', how='left')     
    merge1 = merge1[merge1['Gene Symbol'].notna()] 
    merge1 = merge1.drop_duplicates(subset=['Gene Symbol'])
    intermediate = interpoltools['interpolate'].merge(merge1, on='Gene Symbol', how='left')
    intermediate = intermediate.drop(['Gene Symbol','Human Kinase'], axis = 1)
    intermediate = intermediate.rename(columns = {'Interpolated':'Kinases'})
    merge1 = merge1[['Gene Symbol']
                          +[c for c in merge1 if c not in ['Gene Symbol']]].drop(['Human Kinase'], axis = 1)
    merge1 = merge1.rename(columns = {'Gene Symbol':'Kinases'})
    tot_kin = pd.concat([merge1,intermediate])
    tot_kin = tot_kin.drop_duplicates(subset=['Kinases'])
    kin521 = interpoltools['521 kinases'].merge(tot_kin,on='Kinases',how='left')
    return kin521


# =============================================================================
# calculator for the cMBF
# input listvals is the input data for the new likelihood vector
# noise is intrinsic noise, will depend on if data is expression data,
# dot product scoring, correlation coefficients, etc
# be aware of minimum likelihood cutoff value
# =============================================================================
def cMBF_calc(newdata,noise):
    cMBF = list(range(len(newdata)))
    for i in range(len(newdata)):
        if np.isnan(newdata[i]) == False:
            cMBF[i] = 1-math.exp(-(newdata[i]/noise)*(newdata[i]/noise)/2)
        else:
            cMBF[i] = 0
    return cMBF


# =============================================================================
# single position calculator
# intput: cluster is the cluster (entered as a string eg 'IB') we wish to evaluate and
# dotposdata is the interpolated dot product position-cluster data
# =============================================================================

def selectpos_bayes(cluster,dotposdata,coloc):
    url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/clus_impt_pos.xlsx?raw=true"
    data = rq.get(url).content
    cluspos = pd.read_excel(BytesIO(data), sheet_name = "positioning")
    clusmat = cluspos[cluster]
    int_posselec = []
    int_posselec.append(dotposdata[0]['Kinases'])
    position =  []
    for n in range(len(clusmat)):
        if clusmat[n] == 1:
            int_posselec.append(dotposdata[n][cluster])
            position.append(n)
    
    if len(position) == 3:
        bayes_output = coloc[['Kinases']].copy()  
        first_frame = pd.concat([int_posselec[0].to_frame(),int_posselec[1].to_frame()], axis=1)
        variable = coloc.merge(first_frame,on='Kinases',how='left')
        pos_cMBF = cMBF_calc(variable[cluster],np.mean(variable[cluster][variable[cluster]>=0]))
        prod = pos_cMBF*variable[cluster + '_known']
        bayes_output[cluster + '_' + str(position[0])] = prod/sum(prod)  
    
        sec_frame = pd.concat([int_posselec[0].to_frame(),int_posselec[2].to_frame()], axis=1)
        variable_1 = bayes_output.merge(sec_frame,on='Kinases',how='left')
        pos_cMBF = cMBF_calc(variable_1[cluster],np.mean(variable_1[cluster][variable_1[cluster]>=0]))
        prod = pos_cMBF*variable_1[cluster + '_' + str(position[0])]
        bayes_output[cluster + '_' + str(position[1])] = prod/sum(prod)   
    
        third_frame = pd.concat([int_posselec[0].to_frame(),int_posselec[3].to_frame()], axis=1)
        variable_2 = bayes_output.merge(third_frame,on='Kinases',how='left')
        pos_cMBF = cMBF_calc(variable_2[cluster],np.mean(variable_2[cluster][variable_2[cluster]>=0]))
        prod = pos_cMBF*variable_2[cluster + '_' + str(position[1])]
        bayes_output[cluster + '_' + str(position[2])] = prod/sum(prod)
    elif len(position) == 2:
        bayes_output = coloc[['Kinases']].copy()  
        first_frame = pd.concat([int_posselec[0].to_frame(),int_posselec[1].to_frame()], axis=1)
        variable = coloc.merge(first_frame,on='Kinases',how='left')
        pos_cMBF = cMBF_calc(variable[cluster],np.mean(variable[cluster][variable[cluster]>=0]))
        prod = pos_cMBF*variable[cluster + '_known']
        bayes_output[cluster + '_' + str(position[0])] = prod/sum(prod)  
    
        sec_frame = pd.concat([int_posselec[0].to_frame(),int_posselec[2].to_frame()], axis=1)
        variable_1 = bayes_output.merge(sec_frame,on='Kinases',how='left')
        pos_cMBF = cMBF_calc(variable_1[cluster],np.mean(variable_1[cluster][variable_1[cluster]>=0]))
        prod = pos_cMBF*variable_1[cluster + '_' + str(position[0])]
        bayes_output[cluster + '_' + str(position[1])] = prod/sum(prod)
    elif len(position) == 1:
        bayes_output = coloc[['Kinases']].copy()  
        first_frame = pd.concat([int_posselec[0].to_frame(),int_posselec[1].to_frame()], axis=1)
        variable = coloc.merge(first_frame,on='Kinases',how='left')
        pos_cMBF = cMBF_calc(variable[cluster],np.mean(variable[cluster][variable[cluster]>=0]))
        prod = pos_cMBF*variable[cluster + '_known']
        bayes_output[cluster + '_' + str(position[0])] = prod/sum(prod)
    else:
        bayes_output = coloc[['Kinases']].copy()  
    return bayes_output







        