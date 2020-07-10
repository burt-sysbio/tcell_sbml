# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:26:09 2020

@author: Philipp
"""


import tellurium as te
import pandas as pd
import seaborn as sns
from scipy import constants
sns.set(context = "poster")
from datetime import date
import os
import pickle


def filter_cells(cells, names):
    out = cells[cells.celltype.isin(names)]
    return out

def compute_cell_states(df):
    df["Precursors"] = df.Prec + df.Prec1
    df["Th1_all"] = df.Th1 + df.Th1_2 + df.Th1_mem
    df["Tfh_all"] = df.Tfh + df.Tfh_2 + df.Tfhc + df.Tfhc_2 + df.Tfh_mem
    df["Tr1"] = df.Tr1 + df.Tr1+df.Tr1_2
    df["Total_CD4"] = df.Precursors + df.Th1_all + df.Tfh_all + df.Tr1
    df["Total_eff"] = df.Th1_all + df.Tfh_all + df.Tr1
    df["nonTfh"] = df.Th1_all+df.Tr1
    df["Tfh_chronic"] = df.Tfhc + df.Tfhc_2
    df["Chronic CD4"] = df.Tfh_chronic + df.Tr1

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
    
    cells = filter_cells(df_tidy, ["Th1_all", "Tfh_all", "Tr1", "nonTfh",
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
modelname = "../model_versions/no_cyto_communication.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()

# =============================================================================
# run simulations
# =============================================================================
r = te.loada(antimony_model)

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


#print(te.getODEsFromModel(r))

# =============================================================================
# combine acute and chronic sim
# =============================================================================
arm_df = pd.DataFrame(arm_sim, columns = arm_sim.colnames)
arm_df["Infection"] = "Arm"
cl13_df = pd.DataFrame(cl13_sim, columns = cl13_sim.colnames)
cl13_df["Infection"] = "Cl13"

# tidy data frames sort cells and cytokines
tidy_arm, cells_arm = tidy_sort(arm_df)
tidy_cl13, cells_cl13 = tidy_sort(cl13_df)

tidy_all = pd.concat([tidy_arm, tidy_cl13])

xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0,10,20,30,40,50,60]
g = sns.relplot(data = tidy_all, x = "time", y = "value", col = "celltype",
                col_wrap = 8, hue = "Infection", kind = "line", facet_kws = {"sharey" : False})
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)
#g.savefig(path+today+"/all_species_timecourse.pdf")
