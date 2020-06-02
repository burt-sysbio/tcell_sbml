# -*- coding: utf-8 -*-
"""
Created on Tue May 19 09:59:30 2020

@author: Philipp

load and clean data
"""

import pandas as pd
import numpy as np

time = [9,30,60]

fahey_arm = "fahey_arm.csv"
fahey_cl13 = "fahey_cl13.csv"

df_fahey_arm = pd.read_csv(fahey_arm)
df_fahey_arm["name"] = "Arm"
df_fahey_arm["time"] = time

df_fahey_cl13 = pd.read_csv(fahey_cl13)
df_fahey_cl13["name"] = "Cl13"
df_fahey_cl13["time"] = time

df_fahey = pd.concat([df_fahey_arm, df_fahey_cl13])
df_fahey = df_fahey.rename(columns = {"tfh_y" : "tfh", "nontfh_y" : "non_tfh"})
df_fahey = df_fahey[["tfh", "non_tfh", "name", "time"]]

print(df_fahey)
df_fahey = df_fahey.melt(id_vars = ["name", "time"], var_name = "celltype")


df_fahey.to_csv("fahey_clean.csv")

fahey_rel = "fahey_rel_tfh.csv"
df_fahey_rel = pd.read_csv(fahey_rel)
df_fahey_rel = df_fahey_rel[["arm_rel_tfh_y", "cl13_rel_tfh_y"]]
df_fahey_rel["time"] = time
df_fahey_rel = df_fahey_rel.rename(columns = {"arm_rel_tfh_y" : "Arm", "cl13_rel_tfh_y": "Cl13"})

df_fahey_rel = pd.melt(df_fahey_rel, var_name = "name", id_vars = ["time"], value_name = "rel_tfh")

df_fahey_rel.to_csv("fahey_rel_clean.csv")