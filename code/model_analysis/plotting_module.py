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
sns.set(context = "paper", style = "ticks")


def plot_param_uncertainty(df, pname, save = False): 
    """
    provide data frames generated with run_param_uncertainty from utils
    plots two lineplots on top of each other for armstrong and cl13 
    parameter variation
    """
    
    # split up dataframe for plotting
    arm = df[df.Infection == "Arm"]
    cl13 = df[df.Infection == "Cl13"]
    
    arm_tfh = utils.filter_cells(arm, ["Tfh_all"])
    cl13_tfh = utils.filter_cells(cl13, ["Tfh_all"])
    arm_notfh = utils.filter_cells(arm, ["nonTfh"])
    cl13_notfh = utils.filter_cells(cl13, ["nonTfh"])
    
    # create colorbar
    cmap = "Blues"
    norm = plt.Normalize(arm.val.min(), arm.val.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    # create figure and axes
    fig, axes = plt.subplots(3,2, figsize = (10,10))
    
    ax = axes.flatten()
    
    # tfh relative cell numbers
    sns.lineplot(data = arm_tfh, x = "time", y = "total", hue = "val",
                 palette = cmap, ax = ax[0])
    
    sns.lineplot(data = cl13_tfh, x = "time", y = "total", hue = "val",
                 palette = cmap, ax = ax[1])

    # tfh total cell numbers
    sns.lineplot(data = arm_tfh, x = "time", y = "value", hue = "val",
                 palette = cmap, ax = ax[2])
    
    
    sns.lineplot(data = cl13_tfh, x = "time", y = "value", hue = "val",
                 palette = cmap, ax = ax[3])
    
    # non tfh total cell numbers
    sns.lineplot(data = arm_notfh, x = "time", y = "value", hue = "val",
                 palette = cmap, ax = ax[4])
    
    sns.lineplot(data = cl13_notfh, x = "time", y = "value", hue = "val",
                 palette = cmap, ax = ax[5])
    
   
    reg_scale = [ax[0], ax[1]]
    log_scale = [ax[2], ax[3], ax[4], ax[5]]
    
    ax[0].set_title("Acute Infection")
    ax[1].set_title("Chronic Infection")
    ax[2].set_ylabel("Tfh cells")
    ax[3].set_ylabel("Tfh cells")
    ax[4].set_ylabel("Th1 cells")
    ax[5].set_ylabel("Th1 cells")
            
    for a in reg_scale:
        a.set_ylim(0,100)
        a.set_ylabel("Tfh % of total")
    
    for a in log_scale:
        a.set_yscale("log")
        a.set_ylim(1e0, 1e8)

    for a in ax:
        a.get_legend().remove()
    
    plt.colorbar(sm, ax = [ax[0], ax[1]], label = pname)
    plt.colorbar(sm, ax = [ax[2], ax[3]], label = pname)
    plt.colorbar(sm, ax = [ax[4], ax[5]], label = pname)

    if save == True:
        today = str(date.today())
        path = "../../figures/"
        savepath = path+today
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        fig.savefig(savepath+"/sensitiviy_analysis_"+pname+".pdf")
    # set to origin at the end of experiment
