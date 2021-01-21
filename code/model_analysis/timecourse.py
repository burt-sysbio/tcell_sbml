# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:26:09 2020

@author: Philipp
"""
import tellurium as te
import matplotlib
matplotlib.use("module://backend_interagg")
import pandas as pd
import seaborn as sns
from modules.utils2 import run_pipeline, fit_fun
from datetime import date
import os
import matplotlib.pyplot as plt
sns.set(context = "poster")
# output settings
today = str(date.today())
path = "../../figures/"
if not os.path.exists(path+today):
    os.makedirs(path+today)

# =============================================================================
#load model   
# =============================================================================
modelname = "../models/model_2021.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()

# =============================================================================
# run simulations
# =============================================================================
r = te.loada(antimony_model)
cells, molecules = run_pipeline(r)

# # =============================================================================
# # plot output
# # =============================================================================
xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0,10,20,30,40,50,60]
#
#
# df1 = cells[cells.celltype.isin(["Tfh_all", "nonTfh"])]
# g = sns.relplot(data = df1, x = "time", y = "value", hue = "celltype",
#                 col = "Infection", kind = "line", palette = ["k", "0.5"])
# g.set(yscale = "log", ylabel = ylabel, ylim = (1e2, 1e7), xlim = (0,70),
#       xlabel = xlabel)
# for ax, df in zip(g.axes.flat, data_fahey):
#     sns.scatterplot(data = df, x = "time", y = "value", hue = "celltype",
#                     ax = ax, legend = False, palette = ["0.5", "k"])
#
#     ax.set_ylabel("cells")
#
# # =============================================================================
# # get other output
# # =============================================================================
g = sns.relplot(data = cells, x = "time", y = "value",  hue = "Infection",
                 kind = "line", col = "species", facet_kws = {"sharey": True},
                 col_wrap = 3)
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
plt.show()

g = sns.relplot(data = molecules, x = "time", y = "value",  hue = "Infection",
                 kind = "line", col = "species", facet_kws = {"sharey": False})
g.set(xlabel = xlabel, xticks = xticks)
plt.show()
