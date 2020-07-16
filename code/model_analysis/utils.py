#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 17:13:15 2020

@author: burt
"""
import numpy as np
import pandas as pd
from scipy import constants


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
    df_arm = []
    df_cl13 = []
    
    for val in vals:
        r.resetToOrigin()
        exec("r.%s = %f" % (name, val))
        cells = run_pipeline(r)
        cells_rel = get_rel_cells(cells)
        cells_rel = filter_cells(cells_rel, ["Tfh_all"])
        cells_rel["val"] = val
        
        
        df_arm.append(cells_rel[cells_rel.Infection == "Arm"])
        df_cl13.append(cells_rel[cells_rel.Infection == "Cl13"])
        
    df_arm = pd.concat(df_arm)
    df_cl13 = pd.concat(df_cl13)
    
    sns.lineplot(data = df_arm, x = "time", y = "total", hue = "val",
                 palette = "Blues", ax = ax, legend = False)

    
    sns.lineplot(data = df_cl13, x = "time", y = "total", hue = "val",
                 palette = "Reds", ax = ax, legend = False)
    
    ax.set_xlabel(name)
    ax.set_ylabel("Tfh (% of total)")

    # set to origin at the end of experiment
    r.resetToOrigin()
    


def compute_cell_states(df):
    """
    # takes data frame and computes cell states
    !!!! This might need adjustment depending on the antimony model

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    df["Precursors"] = df.Prec + df.Prec1
    df["Th1_eff"] = df.Th1 + df.Th1_2 + df.Th1_mem
    df["Tfh_eff"] = df.Tfh + df.Tfh_2 + df.Tfh_mem
    df["Tr1_all"] = df.Tr1+df.Tr1_2

    df["nonTfh"] = df.Th1_eff+df.Tr1_all
    df["Tfh_chronic"] = df.Tfhc + df.Tfhc_2
    df["Tfh_all"] = df.Tfh_chronic + df.Tfh_eff
    df["Th_chronic"] = df.Tfh_chronic + df.Tr1_all
    df["Total_CD4"] = df.Precursors + df.Th1_eff + df.Tfh_all + df.Tr1_all
    df["Th_mem"] = df.Th1_mem + df.Tfh_mem
    df["Th_eff"] = df.Th1_eff+ df.Tfh_eff
    return df


def tidy_sort(df):
    # tidy roadrunner sim output
    df = compute_cell_states(df)
    df_tidy = df.melt(id_vars = ["time", "Infection"], var_name = "celltype")
    
    return df_tidy


def get_cytos(df_tidy, keep_cytos = ["IL2", "IL10"]):
    cytos = filter_cells(df_tidy, keep_cytos)
    cytos = cytos.rename(columns = {"celltype" : "cytokine"})  
    cytos["conc_pM"] = get_conc(cytos.value)    
    return cytos


def get_conc(cyto):
    # convert cytokine from molecules to picoMolar
    N = constants.N_A
    # assume lymph node volume is one microlitre
    VOL = 1e-6
    PICOMOL = 1e12

    cyto = cyto*PICOMOL / (VOL*N)
    return cyto


def run_pipeline(r, start = 0, stop = 70, res = 200):
    """
    run both armstrong and cl13 simulation with same parameters
    output cells and cytokines separated

    Parameters
    ----------
    r : TYPE
        DESCRIPTION.
    start : TYPE, optional
        DESCRIPTION. The default is 0.
    stop : TYPE, optional
        DESCRIPTION. The default is 70.
    res : TYPE, optional
        DESCRIPTION. The default is 200.

    Returns
    -------
    cells : TYPE
        DESCRIPTION.
    cytos : TYPE
        DESCRIPTION.

    """
    arm_sim = r.simulate(start, stop, res)
    
    # sum arm and cl13 sim
    r.reset()
    r.deg_TCR = 0.0001
    cl13_sim = r.simulate(start,stop,res)

    # finalze dfs
    arm_df = pd.DataFrame(arm_sim, columns = arm_sim.colnames)
    arm_df["Infection"] = "Arm"
    cl13_df = pd.DataFrame(cl13_sim, columns = cl13_sim.colnames)
    cl13_df["Infection"] = "Cl13"
    
    # tidy data frames sort cells and cytokines
    tidy_arm = tidy_sort(arm_df)
    tidy_cl13 = tidy_sort(cl13_df)
    
    cells = pd.concat([tidy_arm, tidy_cl13])
    cytos = get_cytos(cells)
    return cells, cytos