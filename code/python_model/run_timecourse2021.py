from python_model.parameters import d
from modules.exp_sbml import Sim
from modules.utils2 import run_timecourse
import seaborn as sns
from models.virus_models import vir_model_SIR
import matplotlib.pyplot as plt
import pandas as pd
sns.set(context = "poster", style = "ticks")

# load fit result
fit_name = "20210119_fit"
fit_date = "2021-01-19"
sim = Sim(d, virus_model=vir_model_SIR)
sim.set_fit_params(fit_date, fit_name)

sim.params["SIR_r0"] = 1.2
cells1, molecules1 = run_timecourse(sim, "r1.5", cellnames= "all")
sim.params["SIR_r0"] = 2.5
cells2, molecules2 = run_timecourse(sim, "r2.5", cellnames= "all")
sim.params["SIR_r0"] = 3.0
cells3, molecules3 = run_timecourse(sim, "r10.0", cellnames= "all")

cells = pd.concat([cells1, cells2, cells3])
molecules = pd.concat([molecules1, molecules2, molecules3])

g = sns.relplot(data = cells, x = "time", y = "value", hue = "Infection",
                col = "species", kind = "line", col_wrap= 3)
g.set(yscale = "log", ylim = (1e-1, None))
plt.show()

g = sns.relplot(data = molecules, x = "time", y = "value", hue = "Infection",
                col = "species", kind = "line")
plt.show()