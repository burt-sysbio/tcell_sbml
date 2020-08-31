#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:09:49 2020

@author: burt
sensitivity analysis for parameter annotated model
"""

import tellurium as te
import matplotlib
print(matplotlib.get_backend())
import pandas as pd
import seaborn as sns
sns.set(context = "poster", style = "ticks")
from datetime import date
import os
import matplotlib.pyplot as plt
import utils
import numpy as np

# =============================================================================
# make dir
# =============================================================================
today = str(date.today())
path = "../../figures/"
if not os.path.exists(path+today):
    os.makedirs(path+today)   
# =============================================================================
#load model   
# =============================================================================
#modelname = "../model_versions/20200716_w_cyto_comm.txt"
modelname = "../model_versions/const_precursors.txt"

with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()
    
r = te.loada(antimony_model)
# =============================================================================
# run simulations
# =============================================================================


pnames = ["r_Naive", "death_Prec", "death_Th1", "death_Tfh", "death_Tr1", "death_Tfhc",
          "r_Prec", "n_div", "r_Mem_base", "prolif_TCR_base", "prolif_Th1_base",
          "prolif_Tfh_base", "prolif_Prec_base", "prolif_Tr1_base", "prolif_Tfhc_base",
          "deg_Myc", "EC50_Myc", "deg_TCR", "EC50_TCR"]

celltypes = ["Th1_eff", "Tfh_eff", "Tr1_all", "Tfh_chronic", "Total_CD4"]

reads = utils.sensitivity_analysis2(r, pnames, celltypes)
reads = reads[reads.readout == "Response Size"]
g = sns.catplot(data = reads, x = "pname", y = "log2FC", col = "Infection",
                row = "celltype", kind = "bar", hue = "param_norm")
plt.show()