#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:31:19 2020

@author: burt
"""


import tellurium as te
import pandas as pd
import seaborn as sns

sns.set(context = "poster")
from datetime import date
import os


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
    df["Total_CD4"] = df.Precursors + df.Th1_eff + df.Tfh_all + df.Tr1_all
    df["Th_mem"] = df.Th1_mem + df.Tfh_mem
    df["Th_eff"] = df.Th1_eff+ df.Tfh_eff
    return df


def get_rel_cells(cells):
    #add total cells and compute relative cell fractions
    tot_cells = cells[cells["celltype"] == "Total_CD4"].rename(columns = {"value" : "total"})
    tot_cells = tot_cells[["time", "total", "Infection"]]
    cells = pd.merge(cells, tot_cells, how = "left", on = ["time", "Infection"])
    cells.total = (cells.value / cells.total)*100
    
    return cells


def tidy_sort(df):
    df = compute_cell_states(df)
    df_tidy = df.melt(id_vars = ["time", "Infection"], var_name = "celltype")
    
    cells = filter_cells(df_tidy, ["Th1_eff", "Tfh_all", "Tr1_all", "nonTfh",
                                   "Total_CD4"])    
    return df_tidy, cells

# =============================================================================
# make dir
# =============================================================================
today = str(date.today())
path = "../../figures/"
if not os.path.exists(path+today):
    os.makedirs(path+today)
    
# =============================================================================
#load model   
# =============================================================================
#modelname = "../model_versions/no_cyto_communication.txt"
#with open(modelname, 'r') as myfile:
#    antimony_model = myfile.read()

# =============================================================================
# run simulations
# =============================================================================
#r = te.loada(antimony_model)
#r.integrator = 'gillespie'
#r.integrator.seed = 1234

#arm_sim = r.simulate(0, 70, 200)

# reset variables, keep parameters
#r.reset()
#r.r_IL10_ex = 100
# change TCR signal level so that it does not decrease
# increase rate of IFN-I
#r.deg_TCR = 0.001
#cl13_sim = r.simulate(0,70, 200)

#r.resetToOrigin()




start = 0
stop = 3
res = 30

r = te.loada(
    """
    x' = -k1*x
    x = 1
    k1 = 1
    """)


r.simulate(start, stop, res)
r.plot()
r.reset()

r.gillespie(start, stop, res)
r.plot()