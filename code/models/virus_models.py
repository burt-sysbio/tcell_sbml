from scipy.stats import lognorm
from scipy.interpolate import interp1d
import numpy as np


def sir_parameterization(r0, a, b):
    assert r0>=1
    return 1 / (a*r0**b - 1)


def vir_model_SIR(time, d):
    r0 = d["SIR_r0"]
    SD = sir_parameterization(r0, 1.02, 0.37)
    mean = sir_parameterization(r0, 1.05, 0.18)

    shape, scale = get_lognormdist_params(mean, SD)
    mylognorm = lognorm(s=shape, scale = scale)

    def f(t):
        return mylognorm.pdf(t)
    return f


def get_lognormdist_params(mode, stddev):
    """
    Given the mode and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    p = np.poly1d([1, -1, 0, 0, -(stddev/mode)**2])
    r = p.roots
    sol = r[(r.imag == 0) & (r.real > 0)].real
    shape = np.sqrt(np.log(sol))
    scale = mode * sol
    return shape, scale

def vir_model_const(time, d):
    """
    ag level only depends on vir load, not on time
    """
    time = np.arange(np.min(time), 5 * np.max(time), 0.01)
    s = np.ones_like(time)
    s = s*d["vir_load"]
    f = interp1d(time, s, kind = "zero")
    return f