"""
run simulation and compare with data from fahey et al
"""

import tellurium as te

from datetime import date
import os
import matplotlib
matplotlib.use("TkAgg")
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
# modelname = "../model_versions/20200716_w_cyto_comm.txt"
modelname = "../model_versions/const_precursors.txt"

with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()
# =============================================================================
# run simulations
# =============================================================================
r = te.loada(antimony_model)
use_fit = True
if use_fit:
    with open('fit_result.p', 'rb') as fp:
        fit_result = pickle.load(fp)
        print(fit_result)
    for key, val in fit_result.items():
        r[key]=val
        print(r[key])


cells, cytos = run_pipeline(r)

xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0, 10, 20, 30, 40, 50, 60]

# =============================================================================
# plot fahey data
# =============================================================================
path_data = "../../chronic_infection_model/data_literature/"

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
                col="Infection", kind="line", palette=["k", "0.5"],
                aspect= 0.8)
g.set(yscale="log", ylabel="cells", ylim=(1e4, 3e6), xlim=(0, 70),
      xlabel="time (days)")

g.set_titles("{col_name}")
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
#g.savefig(path+today+"/model_fit_simulation.svg")