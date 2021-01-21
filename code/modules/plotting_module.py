#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 11:37:22 2020

@author: burt
module for plotting
"""
import numpy as np
import os
from datetime import date
import seaborn as sns
import matplotlib.pyplot as plt
import utils
from matplotlib.colors import LogNorm, Normalize

def plot_param_uncertainty(df, pname, log_p = False, save = False): 
    """
    provide data frames generated with run_param_uncertainty from utils
    plots two lineplots on top of each other for armstrong and cl13 
    parameter variation
    """
    def split_df(df):      
    # split up dataframe for plotting
        arm = df[df.Infection == "Arm"]
        cl13 = df[df.Infection == "Cl13"]
        return (arm, cl13)
    
    # these are the cell types that will be plotted
    # Tfh all will be plotted as relative fraction
    cells = ["Tfh_all", "Th1_eff", "Tr1_all", "Tfh_eff", "Tfh_chronic"]
    
    # get all cells and split them up according to arm and cl13
    df_list = [utils.filter_cells(df, [cell]) for cell in cells]
    df_list = [split_df(df) for df in df_list]
    
    # create colorbar mappable either lienar or log scale depending on pparam array
    cmap = "Blues"
    
    if log_p == False:
        norm = Normalize(df.val.min(), df.val.max())
    else:
        norm = LogNorm(df.val.min(), df.val.max())
        
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    # create figure and axes
    fig, axes = plt.subplots(5,2, figsize = (10,14))
    
    # create a panel with arm and cl13 infection and all celltypes as rows
    # create lineplot with color indicating param value
    for df, ax, cell in zip(df_list, axes, cells):
        ylabel = cell

        if cell == "Tfh_all":
            yvar = "total"
        else:
            yvar = "value"

        sns.lineplot(data = df[0], x = "time", y = yvar, hue = "val",
                     palette = cmap, ax = ax[0], hue_norm = norm)
    
        sns.lineplot(data = df[1], x = "time", y = yvar, hue = "val",
                     palette = cmap, ax = ax[1], hue_norm = norm)

        # axes specific changes
        for a in ax:
            if cell == "Tfh_all":
                ylabel = "Tfh % of total"
                a.set_ylim(0,100)
                
            else:
                ylabel = cell
                a.set_yscale("log")
                a.set_ylim(1e0, 1e7)
                
            a.get_legend().remove()
            a.set_ylabel(ylabel)
        
        # add colorbar to each row
        plt.colorbar(sm, ax = [ax[0], ax[1]], label = pname)
            
    ax = axes.flatten()
    ax[0].set_title("Acute Infection")
    ax[1].set_title("Chronic Infection")

    if save == True:
        today = str(date.today())
        path = "../../figures/"
        savepath = path+today
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        fig.savefig(savepath+"/sensitiviy_analysis_"+pname+".pdf")
        plt.close()
    # set to origin at the end of experiment
