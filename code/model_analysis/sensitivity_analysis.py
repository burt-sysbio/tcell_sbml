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


pnames = ["death_Tr1", "death_Tfhc", "r_Prec",
          "prolif_Tr1_base", "prolif_Tfhc_base", "r_Naive", "r_Mem_base"]

celltypes = ["Tfh_all", "nonTfh"]

reads = utils.sensitivity_analysis2(r, pnames, celltypes)
reads = reads[reads.readout == "Response Size"]

g = sns.catplot(data = reads, x = "pname", y = "log2FC", col = "celltype",
                row = "Infection", kind = "bar", hue = "param_norm", palette=["0.2","0.6"], margin_titles= True,
                aspect = 1.5)
g.set(ylabel="Effect Size")
g.set_xticklabels(rotation=90)
plt.show()


for cell in celltypes:
    df = reads[reads.celltype == cell]
    g = sns.catplot(data = df, x = "pname", y = "log2FC", col = "Infection",
                    kind = "bar", hue = "param_norm", palette= ["0.2", "0.6"])
    g.set(ylabel = "Effect Size "+cell)
    g.set_xticklabels(rotation=90)
    #plt.show()
    #g.savefig(path+today+"/sensitivity_" + cell + ".svg")
    #g.savefig(path + today + "/sensitivity_" + cell + ".pdf")


df2 = reads[reads.param_norm == 2.0]
g = sns.catplot(data = df2, x = "pname", y = "log2FC", col = "Infection",
                kind = "bar", hue = "celltype", palette= ["0.2", "0.6"])

#g.set(ylabel="Effect Size", ylim = (-1.6,1.6))
g.set_xticklabels(rotation=90)
#g.savefig(path + today + "/sensitivity_10perc_param_increase.pdf")
#plt.show()