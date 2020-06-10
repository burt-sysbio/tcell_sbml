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

def compute_cell_states(df):
    #df["Precursors"] = df.Precursor + df.Precursor1
    df["Th1_all"] = df.Th1 + df.Th1_noIL2 + df.Th1_mem
    df["Tfh_all"] = df.Tfh + df.Tfh_chronic + df.gcTfh + df.gcTfh_chronic + df.Tfh_mem
    #df["Total_CD4"] = df.Precursors + df.Th1_all + df.Tfh_all + df.Tr1
    #df["Total_eff"] = df.Th1_all + df.Tfh_all + df.Tr1
    df["nonTfh"] = df.Th1_all+df.Tr1
    return df


def transform(sim, data, timepoints = [9, 30, 60]):
    sim = pd.DataFrame(sim, columns = sim.colnames)
    # get simulation data at these time points
    sim = sim[sim.time.isin(timepoints)]
    sim = compute_cell_states(sim)
    # only focus on Tfh and non Tfh
    sim = sim[["Tfh_all", "nonTfh"]]
    sim = sim.melt()
    resid = (data.value.values - sim.value.values) / data.eps.values
    return resid

def tidy(sim, name):
    sim = pd.DataFrame(sim, columns = sim.colnames)
    # get simulation data at these time points
    sim = compute_cell_states(sim)
    # only focus on Tfh and non Tfh
    sim = sim[["time", "Tfh_all", "nonTfh"]]
    sim = sim.melt(id_vars = ["time"])
    sim["Infection"] = name
    return sim

def residual(params, r, data1, data2):
    
    a = params.valuesdict()
    
    for key, val in a.items():
        r[key]=val
    
    tstart = 0
    tstop = 80
    res = 321
    
    sim1 = r.simulate(tstart,tstop,res)
    resid1 = transform(sim1, data1)
    r.reset()
    
    # run model version 2
    # increase rate of IFN-I
    r.deg_TCR = 0.001
    r.r_IFNI = 1
    sim2 = r.simulate(tstart,tstop,res)
    resid2 = transform(sim2, data2)

    resid = np.concatenate((resid1, resid2))
    
    r.resetToOrigin()
    return resid



# =============================================================================
# load literature data and model
# =============================================================================
path_data = "chronic_infection_model/data_literature/"

fahey = "fahey_clean.csv"
df_fahey = pd.read_csv(path_data+fahey)
data_arm = df_fahey[df_fahey.name == "Arm"]
data_arm["eps"] = [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]

data_cl13 = df_fahey[df_fahey.name == "Cl13"]
data_cl13["eps"] = [1e4, 1e4, 1e4, 1e3, 1e3, 1e3]

modelname = "baseprolif_model.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()
r = te.loada(antimony_model)


# =============================================================================
# set parameters
# =============================================================================
params = Parameters()

params.add('deg_Myc', value=0.32, min = 0.2, max = 0.4)
params.add('r_Th1_Tfh_base', value=0.25, min = 0.01, max = 2.0)

params.add('r_Th1_Tr1_base', value=0.25, min = 0.01, max = 2.0)
params.add('r_Tfh_chronic0', value=10.0, min = 0.01, max = 20.0)

params.add("death_chronic", value = 0.01, min = 0, max = 0.24)
params.add('death_Tr1', value=0.07, min = 0.01, max = 0.24)

params.add("prolif_TCR0", value=12.5, min = 1, max = 15)
params.add("prolif_Th10", value=2.3, min = 1, max = 2.8)
params.add("prolif_Tr10", value=1.8, min = 1, max = 2.8)
params.add('prolif_cyto0', value=0.3, min = 0.1, max = 0.5)

params.add('deg_IL2_consumers', value=2000, min = 1000, max = 20000.0)
params.add('deg_IL10_consumers', value=1000, min = 500, max = 10000.0)


#params.add('frequency', value=1.0)


out = minimize(residual, params, args=(r, data_arm, data_cl13))
print(out.message)
out_values = out.params.valuesdict()
print(out_values)
with open('fit_result.p', 'wb') as fit_result:
    pickle.dump(out_values, fit_result, protocol=pickle.HIGHEST_PROTOCOL)


for key, val in out_values.items():
    r[key]=val
# =============================================================================
# run simulations             
# =============================================================================
arm_sim = r.simulate(0, 70, 200)

# reset variables, keep parameters
r.reset()
#r.r_IL10_ex = 100
# change TCR signal level so that it does not decrease
# increase rate of IFN-I
r.deg_TCR = 0.001
r.r_IFNI = 1

cl13_sim = r.simulate(0,70,200)

r.resetToOrigin()


sim1 = tidy(arm_sim, "Arm")
sim2 = tidy(cl13_sim, "Cl13")
cells_tfh = pd.concat([sim1, sim2])

g = sns.relplot(data = cells_tfh, x = "time", y = "value", hue = "variable",
                col = "Infection", kind = "line", palette = ["k", "0.5"])

g.set(yscale = "log", ylabel = "cells", ylim = (1e4, 1e7), xlim = (0,70), 
      xlabel = "time (days)")

data_fahey = [data_arm, data_cl13]

for ax, df in zip(g.axes.flat, data_fahey):
    sns.scatterplot(data = df, x = "time", y = "value", hue = "celltype", 
                    ax = ax, legend = False, palette = ["k", "0.5"])
 
    ax.set_ylabel("cells")
    
g.savefig("tellurium_figs/tcell_modelfit.pdf")