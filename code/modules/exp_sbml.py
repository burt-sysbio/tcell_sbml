import numpy as np
from copy import deepcopy
from scipy.integrate import odeint
from models.python_model import model_2021_py as model2021
import pandas as pd
import pickle

class Sim:
    def __init__(self,
                 params,
                 virus_model,
                 model = model2021,
                 n_states = 15):
        self.model = model
        self.default_params = deepcopy(params)
        self.params = params
        self.n = n_states
        self.virus_model = virus_model

    def init_model(self):
        y0 = np.zeros(self.n)
        y0[0] = self.params["initial_cells"]
        y0[-1] = 1 # myc
        y0[-2] = 1 # myc prec
        return y0

    def simulate(self, start, stop , res):
        time = np.linspace(start, stop, res)
        y0 = self.init_model()
        # initialize virus model for given parameter setting, returns interpolated function
        vir_model = self.virus_model(time, self.params)
        state = odeint(self.model, y0, time, args = (self.params, vir_model))

        colnames = ["Naive", "Prec", "Prec1", "Th1", "Th1_2", "Tfh", "Tfh_2", "Tfhc", "Tfhc_2",
                    "Tr1", "Tr1_2", "Th1_mem", "Tfh_mem", "Myc", "Myc_prec"]

        df = pd.DataFrame(state, columns = colnames)
        df["Ag"] = vir_model(time)
        df["time"] = time
        return df

    def reset(self):
        self.params = deepcopy(self.default_params)

    def set_fit_params(self, fit_date, fit_name):
        # set model parameters to fit
        # load fit result
        fit_dir = "../../output/fit_results/" + fit_date + "/"
        with open(fit_dir + fit_name + '.p', 'rb') as fp:
            fit_result = pickle.load(fp)

        for key, val in fit_result.items():
            self.params[key] = val


