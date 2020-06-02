# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:48:48 2020

@author: Philipp
"""


import tellurium as te
import pandas as pd
import seaborn as sns
from scipy import constants
import numpy as np
import matplotlib.pyplot as plt

sns.set(context = "paper")

def get_rel_cells(cells):
    #add total cells and compute relative cell fractions
    tot_cells = cells[cells["celltype"] == "Total_CD4"].rename(columns = {"value" : "total"})
    tot_cells = tot_cells[["time", "total", "Infection"]]
    cells = pd.merge(cells, tot_cells, how = "left", on = ["time", "Infection"])
    cells.total = (cells.value / cells.total)*100
    
    return cells

def filter_cells(cells, names):
    out = cells[cells.celltype.isin(names)]
    return out

def plot_param_uncertainty(model, startVal, name, ax, num_sims = 20):
    stdDev = 0.5

    # assumes initial parameter estimate as mean and iterates 60% above and below.
    vals = np.linspace((1-stdDev)*startVal, (1+stdDev)*startVal, num_sims)
    for val in vals:
        r.resetToOrigin()
        exec("r.%s = %f" % (name, val))
        cells = run_pipeline(r)
        cells_rel = get_rel_cells(cells)
        cells_rel = filter_cells(cells_rel, ["Tfh_all"])
        
        
        sns.lineplot(data = cells_rel, x = "time", y = "total", hue = "Infection",
                     palette = ["0.2", "tab:red"], ax = ax, legend = False)
        ax.set_xlabel(name)



def compute_cell_states(df):
    df["Precursors"] = df.Precursor + df.Precursor1
    df["Th1_all"] = df.Th1 + df.Th1_noIL2 + df.Th1_mem
    df["Tfh_all"] = df.Tfh + df.Tfh_chronic + df.gcTfh + df.Tfh_mem
    df["Total_CD4"] = df.Precursors + df.Th1_all + df.Tfh_all + df.Tr1
    df["Total_eff"] = df.Th1_all + df.Tfh_all + df.Tr1
    df["nonTfh"] = df.Th1_all+df.Tr1
    return df

def tidy_sort(df):
    
    df = compute_cell_states(df)
    df_tidy = df.melt(id_vars = ["time", "Infection"], var_name = "celltype")
    
    cells = filter_cells(df_tidy, ["Th1_all", "Tfh_all", "Tr1", "nonTfh",
                                   "Total_CD4"])
    
    cytos = filter_cells(df_tidy, ["IL2", "IL10", "IL10_ex"])
    cytos = cytos.rename(columns = {"celltype" : "cytokine"})  
    cytos["conc_pM"] = get_conc(cytos.value)
    
    return df_tidy, cells, cytos



def get_conc(cyto):
    # convert cytokine from molecules to picoMolar
    N = constants.N_A
    # assume lymph node volume is one microlitre
    VOL = 1e-6
    PICOMOL = 1e12

    cyto = cyto*PICOMOL / (VOL*N)
    return cyto

def run_pipeline(r):
    tend = 70
    res = 200
    arm_sim = r.simulate(0, tend, res)
    
    # sum arm and cl13 sim
    r.reset()
    r.p2 = 1000000000
    r.deg_TCR = 0.001
    cl13_sim = r.simulate(0,tend,res)

    # finalze dfs
    arm_df = pd.DataFrame(arm_sim, columns = arm_sim.colnames)
    arm_df["Infection"] = "Arm"
    cl13_df = pd.DataFrame(cl13_sim, columns = cl13_sim.colnames)
    cl13_df["Infection"] = "Cl13"
    
    # tidy data frames sort cells and cytokines
    tidy_arm, cells_arm, cytos_arm = tidy_sort(arm_df)
    tidy_cl13, cells_cl13, cytos_cl13 = tidy_sort(cl13_df)
    
    tidy_all = pd.concat([tidy_arm, tidy_cl13])
    
    cells = pd.concat([cells_arm, cells_cl13])
    return cells


# =============================================================================
# run model    
# =============================================================================
modelname = "baseprolif_model.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()


r = te.loada(antimony_model)

startVals = r.getGlobalParameterValues()
ids = r.getGlobalParameterIds()
pnames = ["fb_stren_pos", "r_Th1_Tfh_base", "deg_IL10",
          "deg_Myc", "EC50_TCR", "prolif_cyto0",
          "r_Th1_Tr1_base", "r_Prec", "r_Th1_noIL2"]

arr = [(a,b) for a,b in zip(startVals, ids) if b in pnames]


fig, axes = plt.subplots(3,3)
for ax, a in zip(axes.flatten(), arr):
    plot_param_uncertainty(r, a[0], a[1], ax)

plt.tight_layout()