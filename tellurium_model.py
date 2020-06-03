# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:46:25 2020

@author: Philipp
"""


import tellurium as te
import pandas as pd
import seaborn as sns
from scipy import constants
sns.set(context = "poster")
from datetime import date
import os

today = str(date.today())
path = "tellurium_figs/"
if not os.path.exists(path+today):
    os.makedirs(path+today)

# =============================================================================
# plotting params
# =============================================================================
xlabel = "time post infection"
ycells = "cells"
ycyto = "conc [pM]"
xticks = [0,10,20,30,40,50,60,70]
# =============================================================================
# functions
# =============================================================================
def filter_cells(cells, names):
    out = cells[cells.celltype.isin(names)]
    return out

def compute_cell_states(df):
    df["Precursors"] = df.Precursor + df.Precursor1
    df["Th1_all"] = df.Th1 + df.Th1_noIL2 + df.Th1_mem
    df["Tfh_all"] = df.Tfh + df.Tfh_chronic + df.gcTfh + df.gcTfh_chronic + df.Tfh_mem
    df["Total_CD4"] = df.Precursors + df.Th1_all + df.Tfh_all + df.Tr1
    df["Total_eff"] = df.Th1_all + df.Tfh_all + df.Tr1
    df["nonTfh"] = df.Th1_all+df.Tr1
    return df

def get_conc(cyto):
    # convert cytokine from molecules to picoMolar
    N = constants.N_A
    # assume lymph node volume is one microlitre
    VOL = 1e-6
    PICOMOL = 1e12

    cyto = cyto*PICOMOL / (VOL*N)
    return cyto

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
    
    cytos = filter_cells(df_tidy, ["IL2", "IL10", "IL10_ex"])
    cytos = cytos.rename(columns = {"celltype" : "cytokine"})  
    cytos["conc_pM"] = get_conc(cytos.value)
    
    return df_tidy, cells, cytos

# =============================================================================
# load literature data
# =============================================================================
path_data = "chronic_infection_model/data_literature/"

fahey = "fahey_clean.csv"
fahey_rel = "fahey_rel_clean.csv"
df_fahey = pd.read_csv(path_data+fahey)
df_fahey_rel = pd.read_csv(path_data+fahey_rel)

# =============================================================================
# load antimony file 
# =============================================================================
modelname = "baseprolif_model.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()


r = te.loada(antimony_model)

# =============================================================================
# write current file to SBML
# =============================================================================


# =============================================================================
# run simulations             
# =============================================================================
arm_sim = r.simulate(0, 70, 200)

# prep model for cl13 sim
r.resetToOrigin()
#r.r_IL10_ex = 100
# change TCR signal level so that it does not decrease
r.deg_TCR = 0.001
# increase rate of IFN-I
r.r_IFNI = 1
cl13_sim = r.simulate(0,70,200)
r.resetToOrigin()


#print(te.getODEsFromModel(r))

# =============================================================================
# finalize data frames
# =============================================================================
arm_df = pd.DataFrame(arm_sim, columns = arm_sim.colnames)
arm_df["Infection"] = "Arm"
cl13_df = pd.DataFrame(cl13_sim, columns = cl13_sim.colnames)
cl13_df["Infection"] = "Cl13"

# tidy data frames sort cells and cytokines
tidy_arm, cells_arm, cytos_arm = tidy_sort(arm_df)
tidy_cl13, cells_cl13, cytos_cl13 = tidy_sort(cl13_df)

tidy_all = pd.concat([tidy_arm, tidy_cl13])
# =============================================================================
# plot everything
# =============================================================================
g = sns.relplot(data = tidy_all, x = "time", y = "value", col = "celltype",
                col_wrap = 8, hue = "Infection", kind = "line", facet_kws = {"sharey" : False})
g.set(yscale = "log", ylim = (0.01, None), xlabel = xlabel, xticks = xticks)


cells = pd.concat([cells_arm, cells_cl13])
cytos = pd.concat([cytos_arm, cytos_cl13])

# =============================================================================
# cl13 vs arm celltypes
# =============================================================================
cells_eff = filter_cells(cells, ["Th1_all", "Tfh_all", "Tr1"])
g = sns.relplot(data = cells_eff, x = "time", y = "value", hue = "Infection",
                col = "celltype", col_wrap = 3, kind = "line")
g.set(yscale = "log", ylabel = "cells", ylim = (1, None), xlabel = xlabel, xticks = xticks)
g.savefig(path+today+"/celltype_kinetics.pdf")

g = sns.relplot(data = cells_eff, x = "time", y = "value", col = "Infection",
                hue = "celltype", kind = "line")
g.set(yscale = "log", ylim = (1, None), xlabel = xlabel, xticks = xticks)

# =============================================================================
# cytokines
# =============================================================================
g = sns.relplot(data = cytos, x = "time", y = "conc_pM", hue = "Infection",
                col = "cytokine", col_wrap = 3, kind = "line",
                facet_kws = {"sharey" : False})


# =============================================================================
# total cells cl13 vs arm
# =============================================================================
cells_tot = cells[cells.celltype == "Total_CD4"]
g = sns.relplot(data = cells_tot, x = "time", y = "value", hue = "Infection",
                kind = "line", aspect = 1.5)

g.set(yscale = "log", ylabel = "total CD4+ cells", ylim = (1e2, None), xlabel = xlabel, xticks = xticks)
g.savefig(path+today+"/total_cells.pdf")
# =============================================================================
# literature figures fahey et al
# =============================================================================
cells_tfh = filter_cells(cells, names = ["Tfh_all", "nonTfh"])

g = sns.relplot(data = cells_tfh, x = "time", y = "value", hue = "celltype",
                col = "Infection", kind = "line", palette = ["k", "0.5"])

g.set(yscale = "log", ylabel = "cells", ylim = (1e3, 1e7), xlim = (0,70), xlabel = xlabel, xticks = xticks)

df_fahey_arm = df_fahey[df_fahey.name == "Arm"]
df_fahey_cl13 = df_fahey[df_fahey.name == "Cl13"]
data_fahey = [df_fahey_arm, df_fahey_cl13]

for ax, df in zip(g.axes.flat, data_fahey):
    sns.scatterplot(data = df, x = "time", y = "value", hue = "celltype", 
                    ax = ax, legend = False, palette = ["k", "0.5"])
 
    ax.set_ylabel("cells")
g.savefig(path+today+"/Tfh_kinetics_fahey.pdf")    
# =============================================================================
# relative cells fahey et al
# =============================================================================
cells_rel = get_rel_cells(cells)
cells_rel = filter_cells(cells_rel, ["Tfh_all"])
g = sns.relplot(data = cells_rel, x = "time", y = "total", hue = "Infection",
                kind = "line", palette = ["0.2", "tab:red"])
g.set(ylabel = "% Tfh/(Th1+Tr1+Tfh)", xlabel = xlabel, xlim = (0, 65), ylim = (0,100),
      xticks = [0,10,20,30,40,50,60])

for ax in g.axes.flat:
    sns.scatterplot(data = df_fahey_rel, x = "time", y = "rel_tfh",
                    hue = "name", ax = ax, palette = ["0.2", "tab:red"], 
                    legend = False)
    ax.set_ylabel("% Tfh/(Th1+Tr1+Tfh)")
    
    
g.savefig(path+today+"/Tfh_rel_fahey.pdf")