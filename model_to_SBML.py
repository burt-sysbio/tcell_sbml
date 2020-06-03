# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 19:35:26 2020

@author: Philipp
"""
import tellurium as te
import os

modelname = "baseprolif_model.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()


r = te.loada(antimony_model)
r.exportToSBML("tcell_sbml.xml")
# =============================================================================
# write current file to SBML
# =============================================================================

