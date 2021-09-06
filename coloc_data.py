# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:43:05 2021

@author: leokt
"""

# =============================================================================
# Data sourced from:
# Yang CR, Tongyoo P, Emamian M, Sandoval PC, Raghuram V, Knepper MA.
# Deep proteomic profiling of vasopressin-sensitive collecting duct 
# cells. I. Virtual Western blots and molecular weight distributions. 
# Am J Physiol Cell Physiol. 2015 Dec 15;309(12):C785-98. 
# doi: 10.1152/ajpcell.00213.2015
# =============================================================================

import pandas as pd
import requests as rq
from io import BytesIO

def coloc_dots(link):
    data = rq.get(link).content
    kinase_coloc = pd.read_excel(BytesIO(data), sheet_name = "kinase_distr")
    cluster_coloc = pd.read_excel(BytesIO(data), sheet_name = "cluster_distr")
    
    cluster_means = cluster_coloc.groupby(["Cluster"]).mean()
    cluster_names = cluster_means.index
    cluster_means = cluster_means.to_numpy()
    
    kinase_names = kinase_coloc['Kinases'].to_numpy()
    kinase_coloc = kinase_coloc.drop(['Kinases'],axis=1)
    kinase_coloc = kinase_coloc.to_numpy()
    
    from bayesian_module import dotscores
    colocdot_df = dotscores(kinase_coloc, cluster_means, kinase_names, cluster_names)
    
    return colocdot_df

