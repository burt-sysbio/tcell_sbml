#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:09:49 2020

@author: burt
sensitivity analysis for parameter annotated model
"""
#import matplotlib
import tellurium as te
import pandas as pd
#matplotlib.use("QT5Agg")
import seaborn as sns
sns.set(context = "poster", style = "ticks")
from datetime import date
import os
import matplotlib.pyplot as plt
import utils
import pickle
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

use_fit = True
if use_fit:
    with open('awesome_fit_result.p', 'rb') as fp:
        fit_result = pickle.load(fp)
        print(fit_result)
    for key, val in fit_result.items():
        r[key]=val
        print(r[key])
# =============================================================================
# run simulations
# =============================================================================
pnames = ["death_Tfh", "death_Tfhc", "prolif_Tr1_base", "prolif_Tfhc_base", "r_Mem_base",
          "pTh1_base", "pTr1_base", "pTfh_base"]

pnames = ["pTh1_base", "pTr1_base", "pTfh_base"]
celltypes = ["Tfh_all", "nonTfh"]

bounds = (0.5,2.0)
reads = utils.sensitivity_analysis2(r, pnames, celltypes, bounds = bounds)
reads = reads[reads.readout == "Response Size"]

reads = reads[reads.param_norm == bounds[1]]
g = sns.catplot(data = reads, x = "pname", y = "log2FC", col = "Infection",
                kind = "bar", hue = "celltype", palette=["0.6","0.2"], margin_titles= True,
                aspect = 1.)
g.set(ylabel="Effect size", ylim = (-1,1), xlabel = "")
g.set_xticklabels(rotation=90)

xlabels = ["pTh1", "pTfh", "pTr1"]
for ax in g.axes.flat:
    ax.set_xticklabels(xlabels)

plt.show()

g.savefig(path+today+"/prob_sensitivity.svg")
g.savefig(path+today+"/prob_sensitivity.pdf")