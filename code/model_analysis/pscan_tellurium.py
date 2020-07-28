# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:48:48 2020

@author: Philipp
"""


import tellurium as te
import os
from datetime import date
from utils import sensitivity_analysis
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(context = "paper", style = "ticks")

today = str(date.today())
path = "../../figures/"
if not os.path.exists(path+today):
    os.makedirs(path+today)


# =============================================================================
# run model    
# =============================================================================
modelname = "../model_versions/20200716_w_cyto_comm.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()


# load model and get parametel values and parameter names
r = te.loada(antimony_model)


# add names for parameter scan
pnames = ["prolif_Th1_base", 
          "prolif_Tfh_base",
          "prolif_TCR_base",
          "prolif_Prec_base",
          "prolif_Tr1_base",
          "prolif_Tfhc_base",
          "deg_Myc",
          "EC50_TCR",
          "death_Prec",
          "death_Th1",
          "death_Tfh",
          "death_Tr1",
          "death_Tfhc",
          "pTh1_base",
          "pTfh_base",
          "pTr1_base",
          "pTfhc_base"
          ]

pnames = ["prolif_Th1_base"]
plt.close("all")
sensitivity_analysis(r, pnames)