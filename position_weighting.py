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


def sugiyama_analysis(link):
    
    import pandas as pd
    import numpy as np
    import requests as rq
    from io import BytesIO
    from bayesian_module import dotscores
    
    #import data
    data = rq.get(link).content
    ks_df = pd.read_excel(BytesIO(data), sheet_name = "Sheet1")
    
    #organize dataframe per kinase
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
    # screen out kinases with low substrate matches (<=20)
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
 
    #import human reference probabilities
    url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/human_PPSP_background.xlsx?raw=true"
    data = rq.get(url).content
    Sbkgd_df = pd.read_excel(BytesIO(data), sheet_name = 'S-center')
    Tbkgd_df = pd.read_excel(BytesIO(data), sheet_name = 'T-center')
    
    Sbkgd_df = Sbkgd_df.drop([20,21])
    Sbkgd_df = Sbkgd_df.drop(['Position:',0,-7,7], axis=1)
    Tbkgd_df = Tbkgd_df.drop([20,21])
    Tbkgd_df = Tbkgd_df.drop(['Position:',0,-7,7], axis=1)
    Sbkgd_array = Sbkgd_df.to_numpy()
    Tbkgd_array = Tbkgd_df.to_numpy()
    #average S and T backgrounds
    STbkgd_array = np.mean([Sbkgd_array,Tbkgd_array],axis=0)/100
    
    #convert observed frequencies to information content scores per kinase
    #IC(a,i)=P(a,i)*log((P(a,i)/Pref(a)),2)
    IC_array = list(range(len(screened_kinases)))
    for b in IC_array:
        IC_array[b] = freq_array_mod[b]*np.log2(freq_array_mod[b]/STbkgd_array)
    
    #extract clusters
    url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/TC_rat_IMCD_clusters.xlsx?raw=true"
    data = rq.get(url).content
    cdf = pd.read_excel(BytesIO(data),None)
    
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
    
    
    #import IMCD rat reference probabilities
    url = "https://github.com/krbyktl/time_course_bayes/blob/master/data_files/TC_rat_IMCD_background.xlsx?raw=true"
    data = rq.get(url).content
    IMCDbkgd_df = pd.read_excel(BytesIO(data))
    IMCDbkgd_array = (IMCDbkgd_df.to_numpy())/100
    
    #create cluster IC matrices
    #IC(a,i)=P(a,i)*log((P(a,i)/Pref(a)),2)
    IC_cluster_array = list(range(len(cluster_array_mod)))
    for b in IC_cluster_array:
        IC_cluster_array[b] = cluster_array_mod[b]*np.log2(cluster_array_mod[b]/IMCDbkgd_array)
    
        
    # %%
    #Creates a dot product scoring between the frequencies of the 0 position in the Sugiyama kinases and the 0 position for the clusters
    STY_freq = list(range(len(freq_ind)))
    for h in range(len(freq_ind)):
        STY_freq[h] = np.transpose(freq_ind[h])[6]
    
    STY_freq_1 = list(range(len(freqclus_ind)))
    for g in range(len(freqclus_ind)):
        STY_freq_1[g] = np.transpose(freqclus_ind[g])[6]
    
    STYdotscores_df = dotscores(STY_freq,STY_freq_1,screened_kinases,clusters_name)
    
    
    # %%
    #Creates dot product scorings between the information content of the positions in the Sugiyama kinases and the positions for the clusters
    
    pos_dotscore_list = list(range(12))
    for j in range(12):
        pos_kin = list(range(len(IC_array)))
        for h in range(len(IC_array)):
            pos_kin[h] = np.transpose(np.nan_to_num(IC_array[h]))[j]
        pos_clus = list(range(len(IC_cluster_array)))
        for g in range(len(IC_cluster_array)):
            pos_clus[g] = np.transpose(np.nan_to_num(IC_cluster_array[g]))[j]
        pos_dotscore_list[j] = dotscores(pos_kin,pos_clus,screened_kinases,clusters_name)
        
    return STYdotscores_df, pos_dotscore_list


