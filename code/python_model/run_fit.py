# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:27:59 2020
@author: Philipp
"""
import numpy as np
from lmfit import minimize, Parameters
import pandas as pd
import os
from modules.utils2 import fit_fun
from datetime import date
import pickle
from python_model.parameters import d
from modules.exp_sbml import Sim
from models.virus_models import vir_model_const

today = str(date.today())
path = "../../output/fit_results/"

if not os.path.exists(path + today):
    os.makedirs(path + today)

# =============================================================================
# get data
# =============================================================================
path_data = "../../data/"
df_fahey = pd.read_csv(path_data + "fahey_data.csv")
data_arm = df_fahey[df_fahey.name == "Arm"]
data_cl13 = df_fahey[df_fahey.name == "Cl13"]

# get model
sim = Sim(d, virus_model= vir_model_const)

# =============================================================================
# set parameters
# =============================================================================
params = Parameters()
params.add('death_tr1', value=0.05, min=0, max=0.2)
params.add('death_tfhc', value=0.01, min=0, max=0.2)
params.add('prolif_tr1', value=2.8, min=2, max=4.0)
params.add('prolif_tfhc', value=4.1, min=3, max=5.0)
params.add("pth1", value=0.06, min=0, max=1.0)
params.add("ptfh", value=0.04, min=0, max=1.0)
params.add("ptr1", value=0.89, min=0, max=1.0)
params.add("ptfhc", expr="1.0-pth1-ptfh-ptr1")
params.add("r_mem", value=0.01, min=0, max=0.2)
params.add("deg_myc", value = 0.32, min =0.28, max = 0.35)
# =============================================================================
# run fitting procedure
# =============================================================================
out = minimize(fit_fun, params, args=(sim, data_arm, data_cl13))
out_values = out.params.valuesdict()
print(out_values)
print(out.message)
print(np.sqrt(out.chisqr))

# store fit result
fit_name = "20210119_fit"
with open(path+today+ "/" + fit_name+'.p', 'wb') as fit_result:
    pickle.dump(out_values, fit_result, protocol=pickle.HIGHEST_PROTOCOL)
