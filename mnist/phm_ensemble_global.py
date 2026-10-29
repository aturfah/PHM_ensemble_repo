import numpy as np
from pathlib import Path
from time import time

from config import SUBSAMP_SIZE, part1_base, r_phm, save_dir_global_ensemble
from helpers import lock_print, load_json, save_json, compute_results
from model_fitting import convert_R_params_object2, robjects, importr, numpy2ri, localconverter, dill
from PHM import PHM

#####################################
### Subagg GMM (mclust All covar) ###
#####################################
def loop_func_ensemble_global(X, truth, idx, lock, pfx, verbose=False, terminate_fitting=False, replicates=50, subsamp_size=SUBSAMP_SIZE):
    if verbose:
        lock_print(lock, idx, "PHM Replicate", idx, flush=True)

    fname = part1_base.format(mthd="PHMAll", rpl=idx)
    ########################################
    #### Load if replicate has been run ####
    ########################################
    if not terminate_fitting:
        if Path(fname).exists():
            if verbose:
                lock_print(lock, idx, "\tLoading from file {}".format(fname), flush=True)
            res = load_json(fname)
            return res

    #################################
    #### Run PHM and cut at K=10 ####
    #################################

    ## Load Densities
    if verbose:
        lock_print(lock, idx, "\tLoading densities...", flush=True)

    prefix = pfx.format(num_pcs=X.shape[1], idx=idx)
    with localconverter(robjects.default_converter + numpy2ri.converter):
        result_part = r_phm.constructPmcParamsSubAggMclust(X, 
                                                           replicates, subsamp_size,
                                                           G=np.array([x for x in range(1, 51)]),
                                                           numCores=1,
                                                           prefix=prefix,
                                                           saveDir=str(save_dir_global_ensemble),
                                                           verbose=False)
    timing_base = "mclustSubagg{cov}{num_pcs}_r{idx}_timing.npy"
    time_density = np.load(
        save_dir_global_ensemble.joinpath(timing_base.format(cov="ALL", num_pcs=X.shape[1], idx=idx))
    )

    if verbose:
        lock_print(lock, idx, "\tLoading complete", flush=True)

    clusters_subagg, _ = convert_R_params_object2(result_part)
    
    if terminate_fitting:
        return clusters_subagg

    # Fit and run PHM
    phm_fname = save_dir_global_ensemble.joinpath(
        "PHM_{}.pkl".format(pfx.format(cov="ALL", num_pcs=X.shape[1], idx=idx))
    )
    phm_time_fname = save_dir_global_ensemble.joinpath(
        "PHM_{}_time.pkl".format(pfx.format(cov="ALL", num_pcs=X.shape[1], idx=idx))
    )
    if phm_time_fname.exists():
        if verbose:
            print("\tLoading PHM from File")
        with phm_fname.open("rb") as sf:
            phm_gauss = dill.load(sf)
        with phm_time_fname.open("rb") as sf:
            phm_gauss_time = dill.load(sf)
    else:
        ## DeltaPmc computation
        if verbose:
            lock_print(lock, idx, "\tComputing deltaPmc", flush=True)
        start = time()
        phm_gauss = PHM(clusters_subagg, "alpha")
        phm_gauss.computeDeltaPmc(mcSamples=10**6, batchSize=10**4)

        ## Run PHM
        if verbose:
            lock_print(lock, idx, "\tRunning PHM", flush=True)
        phm_gauss.runv2()
        phm_gauss_time = time() - start
        with phm_fname.open("wb") as sf:
            dill.dump(phm_gauss, sf)
        with phm_time_fname.open("wb") as sf:
            dill.dump(phm_gauss_time, sf)


    ## Recover K=10 and get labels
    if verbose:
        lock_print(lock, idx, "\tComputing Posterior", flush=True)
    post_mat = phm_gauss.computePosterior(X, 10)
    labels = np.argmax(post_mat, 1).tolist()
    res_repl = compute_results(labels, truth)
    res_repl["time_density"] = [float(time_density)]
    res_repl["time_phm"] = [phm_gauss_time]
    res_repl["time"] = [res_repl["time_phm"][0] + res_repl["time_density"][0]]
    res_repl["n_components"] = [len(clusters_subagg)]
    ## Evaluate
    if verbose:
        lock_print(lock, idx, "\tComplete; Saving {}".format(fname), flush=True)

    try:
        save_json(fname, res_repl)
    except Exception as exc:
        ## Clean up
        Path(fname).unlink()
        raise exc

    return res_repl
