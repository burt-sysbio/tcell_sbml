#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 17:13:15 2020

@author: burt
"""
import numpy as np
import pandas as pd
from scipy import constants
import modules.readout_module as readouts


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


def get_readouts(time, cells):
    """
    get readouts from state array
    """
    peak = readouts.get_peak_height(time, cells)
    area = readouts.get_area(time, cells)
    tau = readouts.get_peaktime(time, cells)
    decay = readouts.get_duration(time, cells)
    
    reads = [peak, area, tau, decay]
    read_names = ["Peak Height", "Response Size", "Peak Time", "Decay"]
    data = {"readout" : read_names, "read_val" : reads}
    reads_df = pd.DataFrame(data = data)
    
    return reads_df


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
        df["Tfh_chr"] = df.Tfhc + df.Tfhc_2
        df["Tfh_all"] = df.Tfh_chr + df.Tfh_eff
        df["Th_chr"] = df.Tfh_chr + df.Tr1_all
        df["Total_CD4"] = df.Precursors + df.Th1_eff + df.Tfh_all + df.Tr1_all
        df["Th_mem"] = df.Th1_mem + df.Tfh_mem
        df["Th_eff"] = df.Th1_eff+ df.Tfh_eff

    else:
        print("no appropriate model name provided")
        
    return df


def get_output(df, name):
    df["Infection"] = name
    df = compute_cell_states(df)
    df_tidy = df.melt(id_vars = ["time", "Infection"], var_name = "species")
    return df_tidy


def get_cytos(df_tidy, keep_cytos = ["IL2", "IL10"]):
    # use again when I use cytos
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


def get_output2(df, cellnames, molnames):
    # should probably check that cellnames and molnames are all in species column.
    cells = df[df.species.isin(cellnames)]
    molecules = df[df.species.isin(molnames)]
    return cells, molecules


def run_timecourse(sim, name, cellnames = "red", molnames = ["Ag", "Myc", "Myc_prec"],
                 start = 0, stop = 70, res = 200):

    sim = sim.simulate(start, stop, res)
    df = get_output(sim, name)
    if cellnames == "all":
        cellnames = ["Precursors", "Th1_eff", "Tfh_eff", "Tr1_all", "nonTfh", "Tfh_chr",
                     "Tfh_all", "Th_chr", "Total_CD4", "Th_mem", "Th_eff"]
    elif cellnames == "red":
        cellnames = ["Precursors", "Th1_eff", "Tfh_eff", "Tfh_chr", "Tr1_all"]

    cells, molecules = get_output2(df, cellnames, molnames)
    return cells, molecules


def run_pipeline(sim, cellnames = "all", molnames = ["Ag", "Myc", "Myc_prec"],
                 start = 0, stop = 70, res = 200):
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
    # run simulation with constant ag for no ag and high ag conditions
    sim.params["vir_load"] = 0
    cells_arm, mols_arm = run_timecourse(sim, name = "Arm", cellnames = cellnames,
                                         molnames = molnames, start = start, stop = stop, res= res)
    
    # sum arm and cl13 sim
    sim.params["vir_load"] = 1
    cells_cl13, mols_cl13 = run_timecourse(sim, name = "Cl13", cellnames = cellnames,
                                         molnames = molnames, start = start, stop = stop, res= res)
    cells = pd.concat([cells_arm, cells_cl13])
    molecules = pd.concat([mols_arm, mols_cl13])
    return cells, molecules


def get_residuals(sim, data, timepoints=[9, 30, 60]):
    """
    formerly called transform
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
    # get simulation data at these time points
    sim = sim[sim.time.isin(timepoints)]
    sim = compute_cell_states(sim)

    # only focus on Tfh and non Tfh
    sim = sim[["Tfh_all", "nonTfh"]]

    # convert to same format as data
    sim = sim.melt()
    resid = (data.value.values - sim.value.values) / data.eps.values
    return resid


def fit_fun(params, r, data1, data2):
    """
    function to be minimized
    run simulation for same parameter set using data1 and data2
    """
    a = params.valuesdict()

    # adjust parammeters
    for key, val in a.items():
        r.params[key] = val

    tstart = 0
    tstop = 80
    res = 321

    # run model arm
    r.params["vir_load"] = 0
    sim1 = r.simulate(tstart, tstop, res)
    resid1 = get_residuals(sim1, data1)

    # run model cl13
    r.params["vir_load"] = 1
    sim2 = r.simulate(tstart, tstop, res)
    resid2 = get_residuals(sim2, data2)

    # array of residuals needs to consist of residuals for Arm and Cl13
    resid = np.concatenate((resid1, resid2))

    r.reset()
    return resid
