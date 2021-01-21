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
from modules.utils2 import run_pipeline
from datetime import date
import os
import pickle
import matplotlib.pyplot as plt

# output settings
sns.set(context = "poster")
today = str(date.today())
path = "../../figures/"
if not os.path.exists(path+today):
    os.makedirs(path+today)

# load data
path_data = "../../data/"
df_fahey = pd.read_csv(path_data + "fahey_data.csv")
data_arm = df_fahey[df_fahey.name == "Arm"]
data_cl13 = df_fahey[df_fahey.name == "Cl13"]

# load fit result
fit_date = today # change this if fit comes from different data than today
fit_dir = "../../output/fit_results/" + fit_date + "/"
fit_name = "2021_fit"
with open(fit_dir + fit_name + '.p', 'rb') as fp:
    fit_result = pickle.load(fp)
    print(fit_result)

# load model
modelname = "../models/model_2021.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()
r = te.loada(antimony_model)

# set model parameters to fit
for key, val in fit_result.items():
    r[key]=val
    print(r[key])

# run simulation
cells, molecules = run_pipeline(r, cellnames= "all")

# =============================================================================
# plot data
# =============================================================================
xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0, 10, 20, 30, 40, 50, 60]

df1 = cells[cells.species.isin(["Tfh_all", "nonTfh"])]
g = sns.relplot(data=df1, x="time", y="value", hue="species",
                col="Infection", kind="line", palette=["k", "0.5"])
g.set(yscale="log", ylabel="cells", ylim=(1e4, 1e7), xlim=(0, 70),
      xlabel="time (days)")

# add fahey data
data_fahey = [data_arm, data_cl13]
for ax, df in zip(g.axes.flat, data_fahey):
    sns.scatterplot(data=df, x="time", y="value", hue="celltype",
                    ax=ax, legend=False, palette=["0.5", "k"])

    ax.set_ylabel("cells")
plt.show()

# g.savefig(path+today+"/modelfit_leastsquares.pdf")

# get other output
g = sns.relplot(data = cells, x = "time", y = "value",  hue = "Infection",
                 kind = "line", col = "species", facet_kws = {"sharey": True}, col_wrap = 3)
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
plt.show()

g = sns.relplot(data = molecules, x = "time", y = "value",  hue = "Infection",
                 kind = "line", col = "species", facet_kws = {"sharey": False})
g.set(xlabel = xlabel, xticks = xticks)
plt.show()