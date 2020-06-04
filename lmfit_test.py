# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:27:59 2020

@author: Philipp
"""


import lmfit

import numpy as np
from lmfit import minimize, Parameters

import tellurium as te
import matplotlib.pyplot as plt

r = te.loada(
    """
    X' = a-b*X
    X = 0
    b = 1
    a = 1
    """)

def residual(params, r, time, data1, data2, eps_data):
    
    a = params.valuesdict()
    
    for key, val in a.items():
        r[key]=val
    
    sim = r.simulate(0,10,20)
    sim_data1 = np.array([sim["X"][np.argmin(np.abs(sim["time"]-t))] for t in time])
    # run model version 2
    r.reset()
    r.X = 0.5
    sim = r.simulate(0,10,20)
    sim_data2 = np.array([sim["X"][np.argmin(np.abs(sim["time"]-t))] for t in time])
    # adjust r for model 2
    
    resid1 = data1-sim_data1
    resid2 = data2-sim_data2
    resid = np.concatenate((resid1, resid2))
    
    r.resetToOrigin()
    return resid


params = Parameters()
params.add('a', value=1.0)
params.add('b', value=0.1)
#params.add('frequency', value=1.0)

time = np.array([1,2,3, 4, 5])
eps_data = np.array([1,1,1,1,1])
data1 = np.array([0.8, 0.9, 1.1, 1.0, 1.2])
data2 = np.array([1.5, 1.7, 1.9, 2.0, 2.0])

out = minimize(residual, params, args=(r, time, data1, data2, eps_data))

out_values = out.params.valuesdict()

for key, val in out_values.items():
    r[key]=val

a = r.simulate(0,5,20)

r.reset()
r.X = 0.5
b = r.simulate(0,5,20)

fig,ax = plt.subplots()
ax.plot(a["time"], a["X"])
ax.plot(b["time"], b["X"])

ax.scatter(time, data1)
ax.scatter(time, data2)