# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 19:40:56 2021

@author: leokt
"""

# =============================================================================
# Data sourced from:
# Sugiyama, N., Imamura, H. & Ishihama, Y. Large-scale Discovery of Substrates 
# of the Human Kinome. Sci Rep 9, 10503 (2019). 
# https://doi.org/10.1038/s41598-019-46385-4
# =============================================================================


import pandas as pd
import numpy as np

#import data and clean up
fileloc = "/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/clean sugiyama kinase substrates.xlsx"
ks_df = pd.read_excel(fileloc, engine='openpyxl')

#divide dataframe per kinase
kinases = ks_df['Kinase'].unique()

AA = ['A','C','D','E','F','G','H','I','J','K','L','M','N','P','Q','R','S','T','V','W','Y']

kin_ind = list(range(len(kinases)))
for i in kin_ind:
    kin_ind[i] = (ks_df[ks_df['Kinase']==kinases[i]])

lengths = list(range(len(kin_ind)))
for i in range(len(lengths)):
    lengths[i] = len(kin_ind[i])
kinase_length_dict = {k:v for k,v in zip(kinases,lengths)}

# %%
# cut out kinases with low peptide matches
screened_kinases = []
for e in range(len(kinase_length_dict)):
    if kinase_length_dict[kinases[e]] > 20:
        screened_kinases.append(kinases[e])

kin_ind_1 = []
for r in range(len(screened_kinases)):
    for w in range(len(kin_ind)):
        if screened_kinases[r] == kin_ind[w]['Kinase'].unique().tolist()[0]:
            kin_ind_1.append(kin_ind[w])

# %%
#function that generates positional frequencies
def make_freq(sort_nam,names):
    freq_matrix = list(range(len(names)))
    for i in range(len(names)):                
        M = 21
        #columns
        N = 13
        init = [ [ 0 for i in range(M) ] for j in range(N) ]        
        ns = (sort_nam[i].groupby(['minus 6']).size()/sort_nam[i].groupby(['minus 6']).size().sum()).to_dict()
        nf = (sort_nam[i].groupby(['minus 5']).size()/sort_nam[i].groupby(['minus 5']).size().sum()).to_dict()
        nfr = (sort_nam[i].groupby(['minus 4']).size()/sort_nam[i].groupby(['minus 4']).size().sum()).to_dict()
        nt = (sort_nam[i].groupby(['minus 3']).size()/sort_nam[i].groupby(['minus 3']).size().sum()).to_dict()
        ntw = (sort_nam[i].groupby(['minus 2']).size()/sort_nam[i].groupby(['minus 2']).size().sum()).to_dict()
        no = (sort_nam[i].groupby(['minus 1']).size()/sort_nam[i].groupby(['minus 1']).size().sum()).to_dict()
        zero = (sort_nam[i].groupby(['zero']).size()/sort_nam[i].groupby(['zero']).size().sum()).to_dict()
        o = (sort_nam[i].groupby(['plus 1']).size()/sort_nam[i].groupby(['plus 1']).size().sum()).to_dict()
        tw = (sort_nam[i].groupby(['plus 2']).size()/sort_nam[i].groupby(['plus 2']).size().sum()).to_dict()
        t = (sort_nam[i].groupby(['plus 3']).size()/sort_nam[i].groupby(['plus 3']).size().sum()).to_dict()
        fr = (sort_nam[i].groupby(['plus 4']).size()/sort_nam[i].groupby(['plus 4']).size().sum()).to_dict()
        f = (sort_nam[i].groupby(['plus 5']).size()/sort_nam[i].groupby(['plus 5']).size().sum()).to_dict()
        s = (sort_nam[i].groupby(['plus 6']).size()/sort_nam[i].groupby(['plus 6']).size().sum()).to_dict()
        freq = [ns,nf,nfr,nt,ntw,no,zero,o,tw,t,fr,f,s]
        for j in freq: 
            for k in range(M):
                if AA[k] not in j:
                    j[AA[k]] = 0
    
    #map frequencies into a 21x13 matrix 
        for l in range(21): 
            for m in range(13):
                init[m][l] = freq[m][AA[l]]      
        print(init)
        freq_matrix[i] =  np.transpose(np.asarray(init))        
    return freq_matrix

freq_ind = make_freq(kin_ind_1,screened_kinases)
    
#exclude "J" and 0 position entry to convert to 20x12 matrix
freq_array_mod = list(range(len(screened_kinases)))
for n in range(len(freq_ind)):
    mod_array = np.delete(freq_ind[n],8,0)
    mod_array = np.delete(mod_array,6,1)
    freq_array_mod[n] = mod_array

#placeholder dictionary
kinase_freq_dict = {k:v for k,v in zip(screened_kinases,freq_array_mod)}

#introduce reference probabilities
fileloc1 = "/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/human PPSP background.xlsx"
Sbkgd_df = pd.read_excel(fileloc1, 'S-center')
Tbkgd_df = pd.read_excel(fileloc1, 'T-center')
#Ybkgd_df = pd.read_excel(fileloc1, 'Y-center')
Sbkgd_df = Sbkgd_df.drop([20,21])
Sbkgd_df = Sbkgd_df.drop(['Position:',0,-7,7], axis=1)
Tbkgd_df = Tbkgd_df.drop([20,21])
Tbkgd_df = Tbkgd_df.drop(['Position:',0,-7,7], axis=1)
Sbkgd_array = Sbkgd_df.to_numpy()
Tbkgd_array = Tbkgd_df.to_numpy()
#average S and T backgrounds
STbkgd_array = np.mean([Sbkgd_array,Tbkgd_array],axis=0)/100

#convert observed frequencies to information content scores [matrix] per kinase
#IC(a,i)=P(a,i)*log((P(a,i)/Pref(a)),2)
IC_array = list(range(len(screened_kinases)))
for b in IC_array:
    IC_array[b] = freq_array_mod[b]*np.log2(freq_array_mod[b]/STbkgd_array)


#extract clusters
fileloc2 = "/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/TC rat IMCD clusters.xlsx"
cdf = pd.read_excel(fileloc2, None)

clusters_name = ['IB','IA1','IA2a','IA2b','IA2c','IIA1a','IIA1b',
            'IIB1a','IIB1b','IIIA','IIIB1','IIIB2','IVA','IVB']

clus_ind = [cdf['I.B'],cdf['I.A.1'],cdf['I.A.2.a'],cdf['I.A.2.b'],
            cdf['I.A.2.c'],cdf['II.A.1.a'],cdf['II.A.1.b'],cdf['II.B.1.a'],
            cdf['II.B.1.b'],cdf['III.A'],cdf['III.B.1'],cdf['III.B.2'],cdf['IV.A'],cdf['IV.B']]

#create frequency matrices
freqclus_ind = make_freq(clus_ind,clusters_name)

#exclude "J" and 0 position entry to convert to 20x12 matrix
cluster_array_mod = list(range(len(freqclus_ind)))
for n in range(len(freqclus_ind)):
    mod_array = np.delete(freqclus_ind[n],8,0)
    mod_array = np.delete(mod_array,6,1)
    cluster_array_mod[n] = mod_array


#import IMCD rat background
fileloc3 = "/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/TC rat IMCD background.xlsx"
IMCDbkgd_df = pd.read_excel(fileloc3)
IMCDbkgd_array = (IMCDbkgd_df.to_numpy())/100

#create cluster IC matrices
#IC(a,i)=P(a,i)*log((P(a,i)/Pref(a)),2)
IC_cluster_array = list(range(len(cluster_array_mod)))
for b in IC_cluster_array:
    IC_cluster_array[b] = cluster_array_mod[b]*np.log2(cluster_array_mod[b]/IMCDbkgd_array)

cluster_IC_dict = {k:v for k,v in zip(clusters_name,IC_cluster_array)}

    
#sum of information content for positions
kinase_sums = []
for v in range(12):
    pos_sums = 0
    for b in range(len(IC_array)):
        pos_sums = pos_sums + sum(abs(np.nan_to_num(np.transpose(IC_array[b])[v])))
    kinase_sums.append(pos_sums)

norm_content = []
for b in kinase_sums:
    norm_content.append(b/np.mean(kinase_sums))
    
kinase_sums.insert(6,0)
positions = ['-6','-5','-4','-3','-2','-1','0','+1','+2','+3','+4','+5','+6']
pos_weight = plt.bar(positions, kinase_sums, color = 'grey')
pos_weight = sns.set_style('white')
pos_weight = sns.set_style('ticks')
pos_weight = sns.despine()
plt.ylabel('Bits', fontsize = 18)
plt.ylim(0,450)
plt.xlabel('Position relative to phosphosite', fontsize = 18)
plt.tight_layout()
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Figures\IC position weight screened N20",
            dpi=600)
    
# %%
#STY screening
STY_freq = list(range(len(freq_ind)))
for h in range(len(freq_ind)):
    STY_freq[h] = np.transpose(freq_ind[h])[6]

STY_freq_1 = list(range(len(freqclus_ind)))
for g in range(len(freqclus_ind)):
    STY_freq_1[g] = np.transpose(freqclus_ind[g])[6]

from bayesmodules import dotscores
STYdotscores_df = dotscores(STY_freq,STY_freq_1,screened_kinases,clusters_name)

fileloc_3 = "/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/STY dot product ranking 355.xlsx"
writer=pd.ExcelWriter(fileloc_3)
STYdotscores_df.to_excel(writer, sheet_name='imported ranking')

writer.save()

# %%
# individual position (IC) screening
pos_dotscore_list = list(range(len(norm_content)))
for j in range(len(norm_content)):
    pos_kin = list(range(len(IC_array)))
    for h in range(len(IC_array)):
        pos_kin[h] = np.transpose(np.nan_to_num(IC_array[h]))[j]
    pos_clus = list(range(len(IC_cluster_array)))
    for g in range(len(IC_cluster_array)):
        pos_clus[g] = np.transpose(np.nan_to_num(IC_cluster_array[g]))[j]
    pos_dotscore_list[j] = dotscores(pos_kin,pos_clus,screened_kinases,clusters_name)
    

# %%
# individual position (freq) screening with weights
pos_dotscorefreq_list = list(range(len(norm_content)))
for j in range(len(norm_content)):
    pos_kin = list(range(len(freq_ind)))
    for h in range(len(freq_ind)):
        pos_kin[h] = np.transpose(np.nan_to_num(freq_ind[h]))[j]
    pos_clus = list(range(len(freqclus_ind)))
    for g in range(len(freqclus_ind)):
        pos_clus[g] = np.transpose(np.nan_to_num(freqclus_ind[g]))[j]
    pos_dotscorefreq_list[j] = dotscores(pos_kin,pos_clus,screened_kinases,clusters_name)
    
# %%
writer=pd.ExcelWriter("/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/position rankings_V3IC.xlsx")
positions = ['-6','-5','-4','-3','-2','-1','+1','+2','+3','+4','+5','+6']
for i, A in enumerate(pos_dotscore_list):
    A.to_excel(writer, sheet_name=positions[i])

writer.save()




