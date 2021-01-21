# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:26:09 2020

@author: Philipp
"""
from python_model.parameters import d
import pandas as pd
import seaborn as sns
from modules.utils2 import run_pipeline
from datetime import date
import os
import pickle
from modules.exp_sbml import Sim
import matplotlib.pyplot as plt
from models.virus_models import vir_model_const
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
fit_name = "20210119_fit"
fit_date = "2021-01-19"

sim = Sim(d, virus_model=vir_model_const)
sim.set_fit_params(fit_date, fit_name)

# run simulation
cells, molecules = run_pipeline(sim, cellnames= "all")

# =============================================================================
# plot data
# =============================================================================
xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0, 10, 20, 30, 40, 50, 60]

df1 = cells[cells.species.isin(["Tfh_all", "nonTfh"])]
g = sns.relplot(data=df1, x="time", y="value", hue="species",
                col="Infection", kind="line", palette=["k", "0.5"])
g.set(yscale="log", ylabel="cells", ylim=(1e3, 1e7), xlim=(0, 70),
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

# get other output
g = sns.relplot(data = molecules, x = "time", y = "value",  hue = "Infection",
                 kind = "line", col = "species", facet_kws = {"sharey": True})
g.set(xlabel = xlabel, xticks = xticks)
plt.show()