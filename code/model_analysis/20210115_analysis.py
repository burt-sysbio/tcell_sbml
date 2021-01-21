"""
run simulation and compare with data from fahey et al
"""
import tellurium as te
import matplotlib
matplotlib.use("module://backend_interagg")
from datetime import date
import os
import matplotlib.pyplot as plt
from utils import run_pipeline, filter_cells
import pandas as pd
import seaborn as sns
import pickle

sns.set(context="poster", style="ticks")
# =============================================================================
# make dir
# =============================================================================
today = str(date.today())
path = "../../figures/"
if not os.path.exists(path + today):
    os.makedirs(path + today)
# =============================================================================
# load model
# =============================================================================
modelname = "../models/model_2021.txt"

with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()
# =============================================================================
# run simulations
# =============================================================================
r = te.loada(antimony_model)

cells, cytos = run_pipeline(r)

xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0, 30, 60]

# =============================================================================
# plot fahey data
# =============================================================================
path_data = "../../data/"

fahey_data = "fahey_cell_numbers.csv"
fahey_error = "fahey_error_bars.csv"

df_fahey = pd.read_csv(path_data + fahey_data)
df_errors = pd.read_csv(path_data + fahey_error)

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
g = sns.relplot(data=df_1, x="time", y="value", hue="celltype",
                row="Infection", kind="line", palette=["k", "0.5"],
                aspect= 0.9, facet_kws= {"sharex" : False})
g.set(yscale="log", ylabel="cells", ylim=(1e4, 3e6), xlim=(0, 70),
      xlabel="time (days)", xticks = xticks)

g.set_titles("{row_name}")

## add literature data
for ax, data, error in zip(g.axes.flat, data_fahey, errors_fahey):
    ax.scatter(data.time, data.iloc[:, 1], c="0.5", s=60)
    ax.scatter(data.time, data.iloc[:, 2], c="k", s=60)

    ax.errorbar(data.time, data.iloc[:, 1], yerr=error.iloc[:, 1], ecolor="0.5",
                fmt="none", capsize=10)
    ax.errorbar(data.time, data.iloc[:, 2], yerr=error.iloc[:, 2], ecolor="k",
                fmt="none", capsize=10)

    ax.set_ylabel("cells")

plt.show()
g.savefig(path+today+"/model_fit_simulation.pdf")
g.savefig(path+today+"/model_fit_simulation.svg")