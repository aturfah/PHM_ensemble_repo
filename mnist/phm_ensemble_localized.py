import numpy as np
from pathlib import Path
from time import time

from config import SUBSAMP_SIZE, part1_base, save_dir_localized_ensemble
from helpers import load_json, save_json, compute_results
from model_fitting import subaggDensityEstimation, leiden_ff, leiden_extract, dill
from PHM import PHM

#####################################
### Subagg GMM (sklearn: k-means) ###
#####################################
def loop_func_ensemble_localized(repl_idx, verbose=False, numCores=1, terminate_fitting=False, replicates=50, subsamp_size=SUBSAMP_SIZE):
    if verbose:
        print(repl_idx, "### REPLICATE {} ###".format(repl_idx))

    ## Load replicate if already ran
    if not terminate_fitting:
        fname = part1_base.format(mthd="PHM_subaggLeiden", rpl=repl_idx)
        print(fname)
        if Path(fname).exists():
            if verbose:
                print(repl_idx, "\tLoading file: {}".format(fname))
            res_repl = load_json(fname)
            res_repl["idx"] = repl_idx
            return res_repl

    ## Otherwise fit    
    if verbose:
        print(repl_idx, "\tDensity Estimation")

    density_time_fname = save_dir_localized_ensemble.joinpath("subaggLeiden_density_timing{}.npy".format(repl_idx))
    start = time()
    res = subaggDensityEstimation(
        X,
        replicates=replicates,
        subSampleSize=subsamp_size,
        fit_func=leiden_ff,
        extract_func=leiden_extract,
        numCores=numCores,
        modelNames=np.asarray(["VVV", "VVI", "VVE"]),
        G=np.asarray([x for x in range(1, 11)]),
        threshold=1e-3,
        verbose=False,
        save_directory=save_dir_localized_ensemble,
        prefix="MNIST_subaggLeiden_R{REP}_{dim}".format(REP=repl_idx, dim=X.shape[1])
    )
    density_time = time() - start
    if density_time_fname.exists():
        density_time = np.load(density_time_fname)
    else:
        np.save(density_time_fname, density_time)

    if terminate_fitting:
        return res
    
    ## Run PHM
    phm_fname = save_dir_localized_ensemble.joinpath("PHM_repl{}.pkl".format(repl_idx))
    phm_time_fname = save_dir_localized_ensemble.joinpath("PHM_repl{}_timing.npy".format(repl_idx))
    if phm_time_fname.exists():
        if verbose:
            print("\tLoading PHM from File")
        phm_time = np.load(phm_time_fname)
        with phm_fname.open("rb") as sf:
            phm_temp = dill.load(sf)
    else:
        if verbose:
            print(repl_idx, "\tComputing deltaPmc")

        start = time()
        phm_temp = PHM(res, "alpha")
        phm_temp.computeDeltaPmc(mcSamples=10**6, batchSize=5 * 10**4)
        if verbose:
            print(repl_idx, "\tRunning PHM")
        phm_temp.runv2()
        phm_time = time() - start

        np.save(phm_time_fname, phm_time)
        with phm_fname.open("wb") as sf:
            dill.dump(phm_temp, sf)

    if verbose:
        print(repl_idx, "\tComputing Posterior Labels")

    post_mat = phm_temp.computePosterior(X, 10)
    labels = np.argmax(post_mat, 1).tolist()

    res_repl = compute_results(labels)
    res_repl["time_density"] = [float(density_time)]
    res_repl["time_phm"] = [float(phm_time)]
    res_repl["time"] = [res_repl["time_density"][0] + res_repl["time_phm"][0]]
    res_repl["n_components"] = [len(res)]
    if verbose:
        print(repl_idx, "\tComplete, saving {}".format(fname))
    try:
        save_json(fname, res_repl)
    except Exception as exc:
        ## Clean up
        Path(fname).unlink()
        raise exc

    res_repl["idx"] = repl_idx
    return res_repl

