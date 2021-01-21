# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:26:09 2020

@author: Philipp
"""


import tellurium as te
import seaborn as sns
import matplotlib
print(matplotlib.get_backend())
sns.set(context = "poster", style = "ticks")
from datetime import date
import os
import matplotlib.pyplot as plt

from utils import run_pipeline, filter_cells
import pandas as pd
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
#modelname = "../models/20200716_w_cyto_comm.txt"
modelname = "../models/const_precursors.txt"

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

fahey_data = "fahey_cell_numbers.csv"
fahey_error = "fahey_error_bars.csv"

df_fahey = pd.read_csv(path_data+fahey_data)
df_errors = pd.read_csv(path_data+fahey_error)

names_arm = ["time", "Tfh_Arm", "NonTfh_Arm"]
names_cl13 = ["time", "Tfh_Cl13", "NonTfh_Cl13"]

data_arm = df_fahey[names_arm]
data_cl13 = df_fahey[names_cl13]
error_arm = df_errors[names_arm]
error_cl13 = df_errors[names_cl13]

data_fahey = [data_arm, data_cl13]
errors_fahey = [error_arm, error_cl13]
# make fahey cell numebrs tidy

# plot Tfh vs non Tfh with literature data
df_1 = filter_cells(cells, ["Tfh_all", "nonTfh"])
g = sns.relplot(data = df_1, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line", palette = ["k", "0.5"])
g.set(yscale = "log", ylabel = "cells", ylim = (1e4, 1e7), xlim = (0,70),
      xlabel = "time (days)")

## add literature data
for ax, data, error in zip(g.axes.flat, data_fahey, errors_fahey):
    ax.scatter(data.time, data.iloc[:,1], c = "0.5", s = 40)
    ax.scatter(data.time, data.iloc[:,2], c = "k", s = 40)
    
    ax.errorbar(data.time, data.iloc[:,1],  yerr= error.iloc[:,1],  ecolor = "0.5", 
                fmt = "none", capsize = 8)
    ax.errorbar(data.time, data.iloc[:,2], yerr = error.iloc[:,2], ecolor = "k", 
                fmt = "none", capsize = 8)
    
    ax.set_ylabel("cells")
    
plt.show()
#g.savefig(path+today+"/model_fit.pdf")

# plot arm vs cl13 per celltype
df_2 = filter_cells(cells, ["Th1_eff", "Tfh_eff", "Tr1_all", "Tfh_chronic",
                               "Precursors", "Total_CD4"])
g = sns.relplot(data = df_2, x = "time", y = "value",  hue = "Infection", 
                kind = "line", col = "celltype", facet_kws = {"sharey": True},
                col_wrap = 3)
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
plt.show()

# plot Tfh vs non Tfh
df_3 = filter_cells(cells, ["Tfh_all", "nonTfh"])
g = sns.relplot(data = df_3, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
plt.show()

# plot chronic vs effector cells
df_4 = filter_cells(cells, ["Th_chronic", "Th_eff"])
g= sns.relplot(data = df_4, x = "time", y = "value", hue = "celltype",
               col = "Infection", kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
plt.show()

# plot memory cells
df_5 = filter_cells(cells, ["Th_mem"])
g= sns.relplot(data = df_5, x = "time", y = "value", hue = "Infection",
               kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
plt.show()

# plot cytokines
g= sns.relplot(data = cytos, x = "time", y = "conc_pM", hue = "Infection",
               col = "cytokine",  kind = "line")
g.set(yscale = "log", ylim = (0.1, None), xlabel = xlabel, xticks = xticks)

plt.show()