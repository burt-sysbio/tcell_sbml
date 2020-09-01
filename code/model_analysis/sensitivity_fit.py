import tellurium as te
import numpy as np
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set(context="poster", style="ticks")
import os
from datetime import date
from utils import compute_cell_states

today = str(date.today())
path = "../../figures/"
if not os.path.exists(path + today):
    os.makedirs(path + today)


def filter_cells(cells, names):
    out = cells[cells.celltype.isin(names)]
    return out


def transform(sim, data, timepoints=[9, 30, 60]):
    """
    objective function computes difference between data and simulation
    Parameters
    ----------
    sim : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    timepoints : TYPE, optional
        DESCRIPTION. The default is [9, 30, 60].

    Returns
    -------
    resid : arr
        array of data-simulation differences (residuals)

    """
    sim = pd.DataFrame(sim, columns=sim.colnames)
    # get simulation data at these time points
    sim = sim[sim.time.isin(timepoints)]
    sim = compute_cell_states(sim)
    # only focus on Tfh and non Tfh
    sim = sim[["Tfh_all", "nonTfh"]]
    # convert to same format as data
    sim = sim.melt()
    resid = (data.value.values - sim.value.values) / data.eps.values
    return resid


def tidy(sim, name):
    """
    convert simulation output to long format for plotting
    """
    sim = pd.DataFrame(sim, columns=sim.colnames)
    # get simulation data at these time points
    sim = compute_cell_states(sim)
    # only focus on Tfh and non Tfh
    sim = sim[["time", "Tfh_all", "nonTfh"]]
    sim = sim.melt(id_vars=["time"])
    sim["Infection"] = name
    return sim


def tidy_sort(df):
    df = compute_cell_states(df)
    df_tidy = df.melt(id_vars=["time", "Infection"], var_name="celltype")

    cells = filter_cells(df_tidy, ["Th1_eff", "Tfh_all", "Tr1_all", "nonTfh",
                                   "Total_CD4"])
    return df_tidy, cells


def residual(params, r, data1, data2):
    """
    function to be minimized
    run simulation for same parameter set using data1 and data2
    """
    r.resetToOrigin()
    a = params.valuesdict()

    # adjust parammeters
    for key, val in a.items():
        r[key] = val

    tstart = 0
    tstop = 80
    res = 321

    # run model arm
    sim1 = r.simulate(tstart, tstop, res)
    resid1 = transform(sim1, data1)
    r.reset()

    # run model cl13
    r.deg_TCR = 0.001
    sim2 = r.simulate(tstart, tstop, res)
    resid2 = transform(sim2, data2)

    # array of residuals needs to consist of residuals for Arm and Cl13
    resid = np.concatenate((resid1, resid2))

    r.resetToOrigin()
    return resid


def fit_sensitivity(r, params, pname, data_arm, data_cl13, bounds):
    """
    take parameter (pname) and vary it, then fit model (provide params dict lmfit object) and return error
    also needs bounds and data for arm and cl13
    """

    assert pname in params.keys()
    myparam = params[pname]
    startVal = myparam.value
    vals = [bounds[0] * startVal, startVal, bounds[1] * startVal]
    err_list = []

    for val in vals:
        # set desired parameter to new value (keep original bounds)
        params[pname].set(value = val, min = myparam.min, max = myparam.max, vary = True)
        print(params)
        # run fit routine, get error (not root calculated by default) and add to list
        fit = minimize(residual, params, args=(r, data_arm, data_cl13))
        err = fit.chisqr
        err = np.sqrt(err)
        err_list.append(err)
        print(pname)
        print(err)
        print(fit.params.valuesdict())

    # set parameters back to original state
    params[pname].set(value=startVal, min=myparam.min, max=myparam.max, vary=True)

    err_arr = np.array(err_list)

    err_norm = err_arr-err_arr[1]
    df = pd.DataFrame({"error": err_list, "error_norm": err_norm})
    df["Param. change"] = [bounds[0], 1.0, bounds[1]]
    df["pname"] = pname

    r.resetToOrigin()
    return df


def run_fit_sensitivity(pnames, r, params, data_arm, data_cl13, bounds = (1.0, 1.0)):
    """
    for each parameter in pnames (list of strings) vary parameter value and fit model afterwords, return error
    """
    df_list = [fit_sensitivity(r, params, pname, data_arm, data_cl13, bounds) for pname in pnames]
    df = pd.concat(df_list)
    return df

# =============================================================================
# load literature data, add uncertainties
# =============================================================================
path_data = "../../chronic_infection_model/data_literature/"

fahey = "fahey_cell_numbers.csv"
df_fahey = pd.read_csv(path_data + fahey)
df_fahey = pd.melt(df_fahey, id_vars=["time"])
df_fahey[["celltype", "name"]] = df_fahey.variable.str.split("_", expand=True)

data_arm = df_fahey[df_fahey.name == "Arm"]
data_arm["eps"] = [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]

data_cl13 = df_fahey[df_fahey.name == "Cl13"]
data_cl13["eps"] = [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]

# load model
modelname = "../model_versions/const_precursors.txt"
with open(modelname, 'r') as myfile:
    antimony_model = myfile.read()

# =============================================================================
# set parameters
# =============================================================================

params = Parameters()
params.add('death_Tr1', value=0.047, min = 0, max = 0.2)
params.add('death_Tfhc', value=0.025, min = 0, max = 0.04)
params.add('prolif_Tr1_base', value=2.9, min = 0.5, max = 6.0)
params.add('prolif_Tfhc_base', value=0.34, min = 0, max = 2.3)
params.add("pTh1_base", value = 0.13, min = 0, max = 1.0)
params.add("pTfh_base", value = 0.08, min = 0, max = 1.0)
params.add("pTr1_base", value = 0.25, min = 0, max = 1.0)
params.add("deg_Myc", value = 0.32, min =0.28, max = 0.35)
params.add("pTfhc_base", expr="1.0-pTh1_base-pTfh_base-pTr1_base")
params.add("r_Mem_base", value = 0.01, min = 0, max = 0.2)

#pnames = ["death_Tr1", "prolif_Tr1_base", "death_Tfhc", "prolif_Tfhc_base",
#          "r_Mem_base", "pTh1_base", "pTfh_base", "pTr1_base"]

pnames = ["death_Tr1"]

r = te.loada(antimony_model)
r.resetToOrigin()
out = minimize(residual, params, args=(r, data_arm, data_cl13))
out_values = out.params.valuesdict()
print(out_values)
print(out.message)
print(np.sqrt(out.chisqr))


#df_err = run_fit_sensitivity(pnames, r, params, data_arm, data_cl13)
#df_err = df_err[df_err["Param. change"] != 1.0]

#g= sns.catplot(data = df_err, x = "pname", y = "error_norm", hue = "Param. change", kind = "bar",
#               palette=["0.2", "0.6"])

#for ax in g.axes.flat:
#    ax.axes.axhline()
#g.set(ylabel="$\Delta$ Error", ylim = (-60, 60))
#g.set_xticklabels(rotation=90)
#plt.show()

#g.savefig(path+today+"/fit_sensitivity.svg")
