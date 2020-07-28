#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 17:13:15 2020

@author: burt
"""
import numpy as np
import pandas as pd
from scipy import constants
from plotting_module import plot_param_uncertainty

def get_rel_cells(cells):
    """
    cells needs to be tidy data frame
    add new column to cells df which gives percentage of total cells
    """
    #add total cells and compute relative cell fractions
    tot_cells = cells[cells["celltype"] == "Total_CD4"].rename(columns = {"value" : "total"})
    tot_cells = tot_cells[["time", "total", "Infection"]]
    cells = pd.merge(cells, tot_cells, how = "left", on = ["time", "Infection"])
    cells.total = (cells.value / cells.total)*100
    
    return cells


def filter_cells(cells, names):
    out = cells[cells.celltype.isin(names)]
    return out


def run_param_uncertainty(r, startVal, name, param_fc, sym, log_p, num_sims = 50):
    """
    

    Parameters
    ----------
    r : roadrunner instance
        antimony T cell model
    startVal : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.
    param_fc : float
        should be greater 1, specifies range of param to be varied as fold-change
    sym : True or "pos" or "neg"
        specify if param should be varied symmetrically with foldchange or only higher or lower than def value 
    log_p : boolean
        set to True if colorbar and sampling should come from log space
    num_sims : TYPE, optional
        DESCRIPTION. The default is 50.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    # assumes initial parameter estimate as mean and iterates 60% above and below.
    if sym == True:
        minval = (1/param_fc)*startVal
        maxval = param_fc*startVal
    elif sym == "pos":
        minval = startVal
        maxval = param_fc*startVal
    else:
        minval = (1/param_fc)*startVal
        maxval = startVal
    
    if log_p == False:
        vals = np.linspace(minval, maxval, num_sims)
    else:
        vals = np.geomspace(minval, maxval, num_sims)
        
    df = []
    
    # for every value in vals arr, change model value and simulate, store result in df_arm and df_cl13
    for val in vals:
        r.resetToOrigin()
        # set model parameter to value
        exec("r.%s = %f" % (name, val))
        
        # simultate model for arm and cl13
        cells, cytos = run_pipeline(r)
        
        cells_rel = get_rel_cells(cells)
        # add current value as val column
        cells_rel["val"] = val
               
        df.append(cells_rel)

    # combine dataframes        
    df= pd.concat(df)

    r.resetToOrigin()

    return df




    
    
def compute_cell_states(df, model_name = "no_cyto_comm_model"):
    """
    # takes data frame and computes cell states
    !!!! This might need adjustment depending on the antimony model

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
        
    model_name : TYPE = String
    provide antimony file name as string because each model has different

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    
    if model_name == "no_cyto_comm_model":
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
    
    else:
        print("no appropriate model name provided")
        
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
    
    r.resetToOrigin()
    
    return cells, cytos


def sensitivity_analysis(r, pnames, param_fc, sym, log_p, save = False):
    """
    run and plot sensitivity analysis for multiple parameters

    Parameters
    ----------
    r : roadrunner instance
        DESCRIPTION.
    pnames : list of parameter names
        DESCRIPTION.

    Returns
    -------
    None

    """
    startVals = r.getGlobalParameterValues()
    ids = r.getGlobalParameterIds()

    # run sensitivity for each parameter name provided
    for val, pname in zip(startVals, ids): 
        if pname in pnames:
            df = run_param_uncertainty(r, val, pname, param_fc, sym, log_p)
            plot_param_uncertainty(df, pname, log_p, save)