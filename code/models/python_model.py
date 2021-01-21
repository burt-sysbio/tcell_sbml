import numpy as np


def pos_fb(cyto, EC50):
    cyto = cyto if cyto > 0 else 0
    out = cyto ** 3 / (cyto ** 3 + EC50 ** 3)
    return out


def neg_fb(cyto, EC50):
    cyto = cyto if cyto > 0 else 0
    out = EC50 ** 3 / (cyto ** 3 + EC50 ** 3)
    return out


def fb_fc(cyto, gamma, EC50):
    cyto = cyto if cyto > 0 else 0
    out = (gamma * cyto ** 3 + EC50 ** 3) / (cyto ** 3 + EC50 ** 3)
    return out


def get_probs(d, ag):
    assert np.abs(d["pth1"] + d["ptfh"] + d["ptr1"] + d["ptfhc"] - 1) < 1e-6
    pTh1 = d["pth1"] * neg_fb(ag, d["EC50_ag"])
    pTfh = d["ptfh"] * neg_fb(ag, d["EC50_ag"])
    pTr1 = d["ptr1"] * pos_fb(ag, d["EC50_ag"])
    pTfhc = d["ptfhc"] * pos_fb(ag, d["EC50_ag"])
    return pTh1, pTfh, pTr1, pTfhc


def norm_probs(probs):
    probs = np.asarray(probs)
    probs = probs / np.sum(probs)
    return probs


def get_eff_prolif(d, myc):
    params = ["prolif_"+ n for n in ["th1", "tfh", "tr1", "tfhc"]]
    prolif_arr = [d[p] * pos_fb(myc, d["EC50_myc"]) for p in params]
    return prolif_arr


def model_2021_py(state, time, d, virus_model):
    naive = state[0]
    prec = state[1]
    prec1 = state[2]
    th1 = state[3]
    th1_2 = state[4]
    tfh = state[5]
    tfh_2 = state[6]
    tfhc = state[7]
    tfhc_2 = state[8]
    tr1 = state[9]
    tr1_2 = state[10]
    th1_mem = state[11]
    tfh_mem = state[12]
    myc = state[13]
    myc_prec = state[14]

    # get probabilities, add feedback by ag and normalize
    ag = virus_model(time)
    p_fb = get_probs(d, ag)
    p_norm = norm_probs(p_fb)
    p_prec = 0.25
    p_norm = (1 - p_prec) * p_norm
    pTh1, pTfh, pTr1, pTfhc = p_norm

    # get proliferation based on myc
    prolif_Th1, prolif_Tfh, prolif_Tr1, prolif_Tfhc = get_eff_prolif(d, myc)
    prolif_Prec = d["prolif_prec"] * pos_fb(myc_prec, d["EC50_myc"])

    d_naive = -d["r_naive"] * naive
    d_prec = 2 * d["r_naive"] * naive + (2 * prolif_Prec * p_prec * d["r_prec"]) * prec1 - (
                d["r_prec"] + d["death_prec"]) * prec

    # precursor differentiation
    d_prec1 = d["r_prec"] * prec - (d["r_prec"] + d["death_prec"]) * prec1

    # effector differentiation
    d_th1 = d["n_div"] * pTh1 * d["r_prec"] * prec1 + prolif_Th1 * (2 * th1_2 - th1) - d["death_th1"] * th1
    d_th1_2 = prolif_Th1 * (th1 - th1_2) - (d["death_th1"] + d["r_mem"]) * th1_2

    d_tfh = d["n_div"] * pTfh * d["r_prec"] * prec1 + prolif_Tfh * (2 * tfh_2 - tfh) - d["death_tfh"] * tfh
    d_tfh_2 = prolif_Tfh * (tfh - tfh_2) - (d["death_tfh"] + d["r_mem"]) * tfh_2

    # chronic Tfh
    d_tfhc = d["n_div"] * pTfhc * d["r_prec"] * prec1 + prolif_Tfhc * (2 * tfhc_2 - tfhc) - d["death_tfhc"] * tfhc
    d_tfhc_2 = prolif_Tfhc * (tfhc - tfhc_2) - d["death_tfhc"] * tfhc_2

    # Tr1 cells
    d_tr1 = d["n_div"] * pTr1 * d["r_prec"] * prec1 + prolif_Tr1 * (2 * tr1_2 - tr1) - d["death_tr1"] * tr1
    d_tr1_2 = prolif_Tr1 * (tr1 - tr1_2) - d["death_tr1"] * tr1_2

    # memory
    d_th1_mem = d["r_mem"] * th1_2 - d["death_mem"]*th1_mem
    d_tfh_mem = d["r_mem"] * tfh_2 - d["death_mem"]*tfh_mem

    # molecules
    deg_myc_prec = d["deg_myc"] * neg_fb(ag, d["EC50_ag"])
    d_myc = -d["deg_myc"] * myc
    d_myc_prec = -deg_myc_prec * myc

    dt_state = [d_naive, d_prec, d_prec1,
                d_th1, d_th1_2, d_tfh, d_tfh_2,
                d_tfhc, d_tfhc_2, d_tr1, d_tr1_2,
                d_th1_mem, d_tfh_mem, d_myc, d_myc_prec]

    return dt_state
