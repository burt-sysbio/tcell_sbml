# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:26:09 2020

@author: Philipp
"""


import tellurium as te
import pandas as pd
import seaborn as sns

sns.set(context = "poster")
from datetime import date
import os
import matplotlib.pyplot as plt
from utils import run_pipeline, filter_cells

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
# =============================================================================
# run simulations
# =============================================================================
r = te.loada(antimony_model)

cells, cytos = run_pipeline(r)


xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0,10,20,30,40,50,60]


# =============================================================================
# plot fahey data
# =============================================================================
plt.close("all")

path_data = "../../chronic_infection_model/data_literature/"

fahey = "fahey_clean.csv"
df_fahey = pd.read_csv(path_data+fahey)
data_arm = df_fahey[df_fahey.name == "Arm"]

data_cl13 = df_fahey[df_fahey.name == "Cl13"]
data_fahey = [data_arm, data_cl13]


# plot Tfh vs non Tfh with literature data
df_1 = filter_cells(cells, ["Tfh_all", "nonTfh"])

g = sns.relplot(data = df_1, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line", palette = ["k", "0.5"])
g.set(yscale = "log", ylabel = "cells", ylim = (1e2, 1e7), xlim = (0,70), 
      xlabel = "time (days)")
for ax, df in zip(g.axes.flat, data_fahey):
    sns.scatterplot(data = df, x = "time", y = "value", hue = "celltype", 
                    ax = ax, legend = False, palette = ["0.5", "k"])
 
    ax.set_ylabel("cells")

#g.savefig(path+today+"/model_fit.pdf")

# plot arm vs cl13 per celltype
df_2 = filter_cells(cells, ["Th1_eff", "Tfh_eff", "Tr1_all", "Tfh_chronic",
                               "Precursors", "Total_CD4"])
g = sns.relplot(data = df_2, x = "time", y = "value",  hue = "Infection", 
                kind = "line", col = "celltype", facet_kws = {"sharey": True},
                col_wrap = 3)
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


# plot Tfh vs non Tfh
df_3 = filter_cells(cells, ["Tfh_all", "nonTfh"])
g = sns.relplot(data = df_3, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


# plot chronic vs effector cells
df_4 = filter_cells(cells, ["Th_chronic", "Th_eff"])
g= sns.relplot(data = df_4, x = "time", y = "value", hue = "celltype",
               col = "Infection", kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


# plot memory cells
df_5 = filter_cells(cells, ["Th_mem"])
g= sns.relplot(data = df_5, x = "time", y = "value", hue = "Infection",
               kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


# plot cytokines
g= sns.relplot(data = cytos, x = "time", y = "conc_pM", hue = "Infection",
               col = "cytokine",  kind = "line")
g.set(yscale = "log", ylim = (0.1, None), xlabel = xlabel, xticks = xticks)

