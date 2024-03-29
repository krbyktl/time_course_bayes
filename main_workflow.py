# -*- coding: utf-8 -*-
"""
Created on Fri May 21 15:51:28 2021

@author: leokt
"""
# =============================================================================
# Steps:

# pre-analyzed:
# mouse IMCD RNA seq (Chen et al 2021)
# IMCD suspension proteome (this paper)
# IMCD microdissected proteome (Limbutara et al 2020)
#
# Screening for S/T kinases (Sugiyama et al 2019)
# Colocalization from fractions (Yang et al 2015)
# Known kinases
# Position ranking (Sugiyama et al 2019) only for selected positions
# =============================================================================

import pandas as pd
import numpy as np
import requests as rq
from io import BytesIO
from bayesian_module import interpolate
from bayesian_module import cMBF_calc
from coloc_data import coloc_dots
from position_weighting import sugiyama_analysis

url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/expression_data.xlsx?raw=true"
data = rq.get(url).content
express_df = pd.read_excel(BytesIO(data), None)
main_vec = express_df['exp_post_vec']

STYdotscores_df, pos_dotscore_list = sugiyama_analysis('https://github.com/krbyktl/time_course_bayes/blob/master/data_files/filtered_kinase_substrates.xlsx?raw=true')
STYdots = interpolate(STYdotscores_df)
main_vec = main_vec.rename(columns = {'Gene Symbol':'Kinases'})


# %%
#STY ranking
main_vec = main_vec.merge(STYdots,on='Kinases',how='left')
clus_names = STYdots.columns.tolist()
clus_names.remove('Kinases')
for i in clus_names: 
    STY_cMBF = cMBF_calc(main_vec[i],1/9)
    prod = STY_cMBF*main_vec['Posterior_Exp']
    main_vec[i + '_STY'] = prod/sum(prod)
main_vec = main_vec.drop(clus_names,axis=1)

# %%
#colocalization ranking
coloc_df = coloc_dots("https://github.com/krbyktl/time_course_bayes/blob/master/data_files/coloc_mapping.xlsx?raw=true")
coloc_df.reset_index(inplace=True)
coloc_df = coloc_df.rename(columns = {'index':'Kinases'})
coloc_vec = main_vec[['Kinases']].copy()
coloc_vec = coloc_vec.merge(coloc_df,on='Kinases',how='left')
for i in clus_names:
    coloc_cMBF = cMBF_calc(coloc_vec[i],np.mean(coloc_vec[i]))
    prod1 = coloc_cMBF*main_vec[i + '_STY']
    coloc_vec[i + '_coloc'] = prod1/sum(prod1)
coloc_vec = coloc_vec.drop(clus_names,axis=1)


# %%
# known kinases
url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/known_kinase_activity.xlsx?raw=true"
data = rq.get(url).content
known_df = pd.read_excel(BytesIO(data), sheet_name = "list")

cluster_direct = ['Increase','Decrease','Decrease','Decrease','Decrease','Increase','Increase','Decrease',
                  'Decrease','Increase','Decrease','Decrease','Increase','Decrease']

clus_names = STYdots.columns.tolist()
clus_names.remove('Kinases')
position_based = coloc_vec
known_based = coloc_vec[['Kinases']].copy()
for a in range(len(clus_names)):
    new_likelihood = list(range(len(position_based['Kinases'])))
    for b in range(len(position_based['Kinases'])):
        if position_based['Kinases'][b] in known_df['Kinases'].values:
            if known_df[known_df['Kinases']==position_based['Kinases'][b]]['Net Effect on Activity'].tolist()[0] == cluster_direct[a]:
                new_likelihood[b] = 0.9
            elif known_df[known_df['Kinases']==position_based['Kinases'][b]]['Net Effect on Activity'].tolist()[0] == 'Regulated':
                new_likelihood[b] = 0.7 
            else:
                new_likelihood[b] = 0.1
        else:
            new_likelihood[b] = 0.5
    known_based[clus_names[a]] = new_likelihood

known_vec = main_vec[['Kinases']].copy()
for c in clus_names:
    prod2 = position_based[c + "_coloc"]*known_based[c]
    new_prior = prod2/sum(prod2)
    known_vec[c+'_known'] = new_prior.tolist()
    
# %%
#position ranking

posit = ['-6','-5','-4','-3','-2','-1','+1','+2','+3','+4','+5','+6']

position_based = main_vec[['Kinases']].copy()
int_pos = list(range(len(pos_dotscore_list)))
for t in range(len(pos_dotscore_list)):
    int_pos[t] = interpolate(pos_dotscore_list[t])

from bayesian_module import selectpos_bayes
from bayesian_module import cMBF_calc
bayes_IB = selectpos_bayes('IB',int_pos,known_vec)

bayes_output = bayes_IB[['Kinases']].copy()  
bayes_output['IB'] = bayes_IB.iloc[:,-1:]

clus_names = STYdots.columns.tolist()
clus_names.remove('Kinases')
rest_clus = clus_names

for v in rest_clus:
    bayes_clus = selectpos_bayes(v,int_pos,known_vec)
    bayes_output[v] = bayes_clus.iloc[:,-1:]
bayes_output['IVA'] = known_vec['IVA_known']
bayes_output['IVB'] = known_vec['IVB_known']



