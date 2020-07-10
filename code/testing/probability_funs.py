# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:48:04 2020

@author: Philipp
"""


import numpy as np
import matplotlib.pyplot as plt

def pos_fb(x, vmax):
    out = x**1/(x**1+1)
    return out*vmax

def neg_fb(x, vmax):
    out = 1/(x**1+1)
    return out*vmax


x = np.arange(0,4,0.1)

th1_0 = 0.4
tfh_0 = 0.4
base = 0.2


th1_arr = pos_fb(x, 0.2)
tfh_arr = neg_fb(x, 0.2)
base_arr = np.ones_like(th1_arr)*base

th1_norm = th1_arr / (th1_arr+tfh_arr+base_arr)
tfh_norm = tfh_arr / (th1_arr+tfh_arr+base_arr)
base_norm = base_arr / (th1_arr+tfh_arr+base_arr)

tot = th1_norm+tfh_norm+base_norm

fig, ax = plt.subplots()
ax.plot(x, th1_norm)
ax.plot(x, tfh_norm)
ax.plot(x, base_norm)
ax.plot(x, tot)