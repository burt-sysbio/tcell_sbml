# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:26:09 2020

@author: Philipp
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

cells_tfh = pd.concat([tidy_arm, tidy_cl13])

xlabel = "time post infection (d)"
ylabel = "cells"
xticks = [0,10,20,30,40,50,60]


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

# =============================================================================
# plot fahey data
# =============================================================================
path_data = "../../chronic_infection_model/data_literature/"

fahey = "fahey_clean.csv"
df_fahey = pd.read_csv(path_data+fahey)
data_arm = df_fahey[df_fahey.name == "Arm"]

data_cl13 = df_fahey[df_fahey.name == "Cl13"]
data_fahey = [data_arm, data_cl13]

df_1 = filter_cells(cells_tfh, ["Tfh_all", "nonTfh"])

g = sns.relplot(data = df_1, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line", palette = ["k", "0.5"])
g.set(yscale = "log", ylabel = "cells", ylim = (1e2, 1e7), xlim = (0,70), 
      xlabel = "time (days)")
for ax, df in zip(g.axes.flat, data_fahey):
    sns.scatterplot(data = df, x = "time", y = "value", hue = "celltype", 
                    ax = ax, legend = False, palette = ["0.5", "k"])
 
    ax.set_ylabel("cells")

g.savefig(path+today+"/model_fit.pdf")
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

