# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:27:59 2020

@author: Philipp
"""


import lmfit

import numpy as np
from lmfit import minimize, Parameters

import tellurium as te
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import seaborn as sns
import os
from datetime import date

today = str(date.today())
path = "../../figures/"
if not os.path.exists(path+today):
    os.makedirs(path+today)
    
    
def filter_cells(cells, names):
    out = cells[cells.celltype.isin(names)]
    return out


def compute_cell_states(df):
    df["Precursors"] = df.Prec + df.Prec1
    df["Th1_eff"] = df.Th1 + df.Th1_2 + df.Th1_mem
    df["Tfh_eff"] = df.Tfh + df.Tfh_2 + df.Tfh_mem
    df["Tr1_all"] = df.Tr1+df.Tr1_2

    df["nonTfh"] = df.Th1_eff+df.Tr1_all
    df["Tfh_chronic"] = df.Tfhc + df.Tfhc_2
    df["Tfh_all"] = df.Tfh_chronic + df.Tfh_eff
    df["Th_chronic"] = df.Tfh_chronic + df.Tr1_all
    df["Th_eff"] = df.Th1_eff + df.Tfh_eff
    df["Th_mem"] = df.Th1_mem+df.Tfh_mem
    df["Total_CD4"] = df.Precursors + df.Th1_eff + df.Tfh_all + df.Tr1_all
    return df


def transform(sim, data, timepoints = [9, 30, 60]):
    """
    objective function computes difference between data and simulation
    Parameters
    ----------
    sim : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    timepoints : TYPE, optional
        DESCRIPTION. The default is [9, 30, 60].

    Returns
    -------
    resid : arr
        array of data-simulation differences (residuals)

    """
    sim = pd.DataFrame(sim, columns = sim.colnames)
    # get simulation data at these time points
    sim = sim[sim.time.isin(timepoints)]
    sim = compute_cell_states(sim)
    # only focus on Tfh and non Tfh
    sim = sim[["Tfh_all", "nonTfh"]]
    # convert to same format as data
    sim = sim.melt()
    resid = (data.value.values - sim.value.values) / data.eps.values
    return resid


def tidy(sim, name):
    """
    convert simulation output to long format for plotting
    """
    sim = pd.DataFrame(sim, columns = sim.colnames)
    # get simulation data at these time points
    sim = compute_cell_states(sim)
    # only focus on Tfh and non Tfh
    sim = sim[["time", "Tfh_all", "nonTfh"]]
    sim = sim.melt(id_vars = ["time"])
    sim["Infection"] = name
    return sim


def tidy_sort(df):
    df = compute_cell_states(df)
    df_tidy = df.melt(id_vars = ["time", "Infection"], var_name = "celltype")
    
    cells = filter_cells(df_tidy, ["Th1_eff", "Tfh_all", "Tr1_all", "nonTfh",
                                   "Total_CD4"])    
    return df_tidy, cells

def residual(params, r, data1, data2):
    """
    function to be minimized
    run simulation for same parameter set using data1 and data2
    """
    a = params.valuesdict()
    
    # adjust parammeters
    for key, val in a.items():
        r[key]=val
    
    tstart = 0
    tstop = 80
    res = 321
    
    # run model arm
    sim1 = r.simulate(tstart,tstop,res)
    resid1 = transform(sim1, data1)
    r.reset()
    
    # run model cl13
    r.deg_TCR = 0.001
    sim2 = r.simulate(tstart,tstop,res)
    resid2 = transform(sim2, data2)

    # array of residuals needs to consist of residuals for Arm and Cl13
    resid = np.concatenate((resid1, resid2))
    
    r.resetToOrigin()
    return resid


# =============================================================================
# load literature data and model
# =============================================================================
path_data = "../../chronic_infection_model/data_literature/"

fahey = "fahey_clean.csv"
df_fahey = pd.read_csv(path_data+fahey)
data_arm = df_fahey[df_fahey.name == "Arm"]
data_arm["eps"] = [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]

data_cl13 = df_fahey[df_fahey.name == "Cl13"]
data_cl13["eps"] = [1e4, 1e4, 1e4, 1e3, 1e3, 1e3]

modelname = "../model_versions/no_cyto_communication.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()
    
r = te.loada(antimony_model)

# =============================================================================
# set parameters
# =============================================================================
params = Parameters()

params.add('death_Th1', value=0.24, min = 0.21, max = 0.27)
params.add('death_Tfh', value=0.24, min = 0.21, max = 0.27)
params.add('death_Tr1', value=0.04, min = 0, max = 0.1)
params.add('death_Tfhc', value=0, min = 0, max = 0.1)
params.add('prolif_Th1_base', value=4.5, min = 4.2, max = 4.7)
params.add('prolif_Tfh_base', value=4.5, min = 4.2, max = 4.7)
params.add('prolif_Tr1_base', value=3.0, min = 2.8, max = 3.2)
params.add('prolif_Tfhc_base', value=2.1, min = 1.9, max = 2.3)
params.add('r_Mem_base', value=0.01, min = 0.005, max = 0.015)


# =============================================================================
# run fitting procedure
# =============================================================================
r.resetToOrigin()
out = minimize(residual, params, args=(r, data_arm, data_cl13))
out_values = out.params.valuesdict()
print(out_values)
print(out.message)
# store fit result
with open('fit_result.p', 'wb') as fit_result:
    pickle.dump(out_values, fit_result, protocol=pickle.HIGHEST_PROTOCOL)



# =============================================================================
# run simulations             
# =============================================================================
# adjust simulation params to fit values
r.resetToOrigin()
for key, val in out_values.items():
    r[key]=val

arm_sim = r.simulate(0, 70, 200)
# reset variables, keep parameters
r.reset()
# change TCR signal level so that it does not decrease
r.deg_TCR = 0.001
cl13_sim = r.simulate(0,70,200)

# =============================================================================
# convert data
# =============================================================================
arm_df = pd.DataFrame(arm_sim, columns = arm_sim.colnames)
arm_df["Infection"] = "Arm"
cl13_df = pd.DataFrame(cl13_sim, columns = cl13_sim.colnames)
cl13_df["Infection"] = "Cl13"

# tidy data frames sort cells and cytokines
tidy_arm, cells_arm = tidy_sort(arm_df)
tidy_cl13, cells_cl13 = tidy_sort(cl13_df)

cells_tfh = pd.concat([tidy_arm, tidy_cl13])
# =============================================================================
# plot data
# =============================================================================
xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0,10,20,30,40,50,60]


df_1 = filter_cells(cells_tfh, ["Tfh_all", "nonTfh"])

g = sns.relplot(data = df_1, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line", palette = ["k", "0.5"])

g.set(yscale = "log", ylabel = "cells", ylim = (1e2, 1e7), xlim = (0,70), 
      xlabel = "time (days)")

data_fahey = [data_arm, data_cl13]

for ax, df in zip(g.axes.flat, data_fahey):
    sns.scatterplot(data = df, x = "time", y = "value", hue = "celltype", 
                    ax = ax, legend = False, palette = ["0.5", "k"])
 
    ax.set_ylabel("cells")
    
g.savefig(path+today+"/modelfit_leastsquares.pdf")

# =============================================================================
# get other output
# =============================================================================
df_2 = filter_cells(cells_tfh, ["Th1_eff", "Tfh_eff", "Tr1_all", "Tfh_chronic",
                               "Precursors", "Total_CD4"])
g = sns.relplot(data = df_2, x = "time", y = "value",  hue = "Infection", 
                kind = "line", col = "celltype", facet_kws = {"sharey": True},
                col_wrap = 3)
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


df_3 = filter_cells(cells_tfh, ["Tfh_all", "nonTfh"])
g = sns.relplot(data = df_3, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


df_4 = filter_cells(cells_tfh, ["Th_chronic", "Th_eff"])
g= sns.relplot(data = df_4, x = "time", y = "value", hue = "celltype",
               col = "Infection", kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


df_5 = filter_cells(cells_tfh, ["Th_mem"])
g= sns.relplot(data = df_5, x = "time", y = "value", hue = "Infection",
               kind = "line")
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
