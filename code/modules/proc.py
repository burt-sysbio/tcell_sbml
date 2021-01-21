def readout_sensitivity(r, startVal, pname, celltypes, bounds):
    """

    Parameters
    ----------
    r : antimony model
    startVal : float
        default parameter value
    pname : string
        parameter name
    celltypes : list
        list of celltype names. e.g. ["Th1_all", "Tfh_all"] as given by compute cell states fun
    bounds : tuple
        increase/decrease default param by this foldchange. The default is (0.9,1.1).

    Returns
    -------
    df : normalized data frame with readouts for different cell types and arm+cl13 infection

    """
    vals = [bounds[0] * startVal, startVal, bounds[1] * startVal]
    reads_list = []

    for val in vals:
        r.resetToOrigin()
        exec("r.%s = %f" % (pname, val))
        start = 0
        stop = 70
        res = 200
        arm_sim = r.simulate(start, stop, res)

        # sum arm and cl13 sim
        r.reset()
        r.deg_TCR = 0.0001
        cl13_sim = r.simulate(start, stop, res)

        # finalze dfs
        arm_df = pd.DataFrame(arm_sim, columns=arm_sim.colnames)
        cl13_df = pd.DataFrame(cl13_sim, columns=cl13_sim.colnames)

        arm_df = compute_cell_states(arm_df)
        cl13_df = compute_cell_states(cl13_df)

        reads = []
        labels = ["Arm", "Cl13"]
        for df, label in zip([arm_df, cl13_df], labels):
            for celltype in celltypes:
                df_reads = get_readouts(df.time, df[celltype])
                df_reads["celltype"] = celltype
                df_reads["Infection"] = label
                reads.append(df_reads)

        reads = pd.concat(reads)
        reads["param_val"] = val

        reads_list.append(reads)

    reads = pd.concat(reads_list)
    reads["pname"] = pname

    df = norm_sens_ana(reads, startVal)
    return df


def run_param_uncertainty(r, startVal, name, param_fc, sym, log_p, num_sims=50):
    """
    Parameters
    ----------
    r : roadrunner instance
        antimony T cell model
    startVal : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.
    param_fc : float
        should be greater 1, specifies range of param to be varied as fold-change
    sym : True or "pos" or "neg"
        specify if param should be varied symmetrically with foldchange or only higher or lower than def value
    log_p : boolean
        set to True if colorbar and sampling should come from log space
    num_sims : TYPE, optional
        DESCRIPTION. The default is 50.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    # assumes initial parameter estimate as mean and iterates 60% above and below.
    if sym == True:
        minval = (1 / param_fc) * startVal
        maxval = param_fc * startVal
    elif sym == "pos":
        minval = startVal
        maxval = param_fc * startVal
    else:
        minval = (1 / param_fc) * startVal
        maxval = startVal

    if log_p == False:
        vals = np.linspace(minval, maxval, num_sims)
    else:
        vals = np.geomspace(minval, maxval, num_sims)

    df = []

    # for every value in vals arr, change model value and simulate, store result in df_arm and df_cl13
    for val in vals:
        r.resetToOrigin()
        # set model parameter to value
        exec("r.%s = %f" % (name, val))

        # simultate model for arm and cl13
        cells, cytos = run_pipeline(r)

        cells_rel = get_rel_cells(cells)
        # add current value as val column
        cells_rel["val"] = val

        df.append(cells_rel)

    # combine dataframes
    df = pd.concat(df)

    r.resetToOrigin()

    return df


def sensitivity_analysis(r, pnames, param_fc, sym, log_p, save=False):
    """
    run and plot sensitivity analysis for multiple parameters

    Parameters
    ----------
    r : roadrunner instance
        DESCRIPTION.
    pnames : list of parameter names
        DESCRIPTION.

    Returns
    -------
    None

    """
    startVals = r.getGlobalParameterValues()
    ids = r.getGlobalParameterIds()

    # run sensitivity for each parameter name provided
    for val, pname in zip(startVals, ids):
        if pname in pnames:
            df = run_param_uncertainty(r, val, pname, param_fc, sym, log_p)
            plot_param_uncertainty(df, pname, log_p, save)


def sensitivity_analysis2(r, pnames, celltypes, bounds):
    """
    run and plot sensitivity analysis for multiple parameters

    Parameters
    ----------
    r : roadrunner instance
        DESCRIPTION.
    pnames : list of parameter names
        DESCRIPTION.

    Returns
    -------
    None
    """
    startVals = r.getGlobalParameterValues()
    ids = r.getGlobalParameterIds()

    # run sensitivity for each parameter name provided
    df_list = []
    for val, pname in zip(startVals, ids):
        if pname in pnames:
            df = readout_sensitivity(r, val, pname, celltypes, bounds)
            df_list.append(df)

    df = pd.concat(df_list)
    # remove the normalized parameter value with FC=0
    df = df[df.param_norm != 1.0]
    return df


def norm_sens_ana(df, norm_val):
    """
    take data frame with readouts from fun readout_sensitivity
    and normalize to middle of 3 values
    """
    df2 = df[df.param_val == norm_val]
    df2 = df2.rename(columns={"read_val": "read_norm"})
    df2 = df2.drop(columns=["param_val"])
    df = df.merge(df2, on=['readout', "celltype", 'pname', "Infection"], how='left')

    df["param_norm"] = norm_val
    df["param_norm"] = np.round(df.param_val / df.param_norm, 1)

    # compute log2FC
    logseries = df["read_val"] / df["read_norm"]
    logseries = logseries.astype(float)

    df["log2FC"] = np.log2(logseries)
    df = df.drop(columns=["read_norm"])
