#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:09:49 2020

@author: burt
sensitivity analysis for parameter annotated model
"""

import tellurium as te
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
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


pnames = ["death_Tr1", "death_Tfhc", "death_Th1", "death_Tfh", "r_Prec", "prolif_Th1_base",
          "prolif_Tr1_base", "prolif_Tfh_base",
          "prolif_Tfhc_base", "prolif_Prec_base"]

celltypes = ["Tfh_all", "nonTfh"]

reads = utils.sensitivity_analysis2(r, pnames, celltypes)
reads = reads[reads.readout == "Response Size"]

for cell in celltypes:
    df = reads[reads.celltype == cell]
    g = sns.catplot(data = df, x = "pname", y = "log2FC", col = "Infection",
                    kind = "bar", hue = "param_norm", palette= ["0.2", "0.6"])
    g.set(ylabel = "Effect Size "+cell)
    g.set_xticklabels(rotation=90)
    plt.show()
    g.savefig(path+today+"/sensitivity_" + cell + ".svg")
    g.savefig(path + today + "/sensitivity_" + cell + ".pdf")



df2 = reads[reads.param_norm == 1.1]
g = sns.catplot(data = df2, x = "pname", y = "log2FC", col = "Infection",
                kind = "bar", hue = "celltype", palette= ["0.2", "0.6"])

g.set(ylabel="Effect Size", ylim = (-1.6,1.6))
g.set_xticklabels(rotation=90)
g.savefig(path + today + "/sensitivity_10perc_param_increase.pdf")
plt.show()