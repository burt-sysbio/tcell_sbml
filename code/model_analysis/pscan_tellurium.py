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
import os
from datetime import date
import matplotlib

sns.set(context = "paper", style = "ticks")

today = str(date.today())
path = "tellurium_figs/"
if not os.path.exists(path+today):
    os.makedirs(path+today)


# =============================================================================
# run model    
# =============================================================================
modelname = "baseprolif_model.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()


r = te.loada(antimony_model)

startVals = r.getGlobalParameterValues()
ids = r.getGlobalParameterIds()
pnames = ["fb_stren_pos", "r_Th1_Tfh_base", "deg_IL10_consumers",
          "prolif_Tr10", "prolif_cyto0", "prolif_Tfh0",
          "r_Th1_Tr1_base", "r_Prec", "r_Th1_noIL2",
          "prolif_Th10", "deg_IL2_consumers", "r_Tfh_chronic0"]


arr = [(a,b) for a,b in zip(startVals, ids) if b in pnames]


fig, axes = plt.subplots(4,3)
for ax, a in zip(axes.flatten(), arr):
    plot_param_uncertainty(r, a[0], a[1], ax)

plt.tight_layout()

fig.savefig(path+today+"/pscan.pdf")



test = plot_param_uncertainty2(r,  0.25 , "r_Prec")



