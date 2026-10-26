import numpy as np
from pathlib import Path
from time import time

from config import part1_base, save_dir_phm_single
from helpers import load_json, save_json, compute_results
from model_fitting import dill, robjects, convert_R_params_object2
from PHM import PHM

def loop_func_phm_single(X, truth, random_state, pfx):
    ## Check if we've already run
    fname = part1_base.format(mthd="PHM_gmmHier", rpl=random_state)
    if Path(fname).exists():
        out = load_json(fname)
        return out

    ## Load/prepare the density
    dens_fname = "{base}/{pfx}{dim}_G50_r{rs}.RDS".format(base=save_dir_phm_single, 
                                                   pfx=pfx, 
                                                   dim=X.shape[1], 
                                                   rs=random_state)
    dens_time_fname = dens_fname.replace(".RDS", "_timing.npy")
    if not Path(dens_fname).exists():
        raise RuntimeError("File Not Found: {}".format(dens_fname))

    if not Path(dens_time_fname).exists():
        raise RuntimeError("File Not Found: {}".format(dens_time_fname))

    clusters, _ = convert_R_params_object2(robjects.r["readRDS"](dens_fname))
    elapsed_density = float(np.load(dens_time_fname))

    ## Run PHM
    phm_fname = Path(dens_fname.replace(".RDS", "_PHM.pkl"))
    phm_time_fname = Path(dens_fname.replace(".RDS", "_PHM_time.pkl"))
    if phm_time_fname.exists():
        with phm_fname.open("rb") as sf:
            phm_mclust = dill.load(sf)
        with phm_time_fname.open("rb") as sf:
            elapsed_phm = dill.load(sf)
    else:
        start = time()
        phm_mclust = PHM(clusters, "alpha")
        phm_mclust.computeDeltaPmc(mcSamples=10**6, batchSize=10**4)
        phm_mclust.runv2()
        elapsed_phm = time() - start
        with phm_fname.open("wb") as sf:
            dill.dump(phm_mclust, sf)

        with phm_time_fname.open("wb") as sf:
            dill.dump(elapsed_phm, sf)

    ## Cut tree at K=10 and finish
    post_mat = phm_mclust.computePosterior(X, 10)
    labels = np.argmax(post_mat, 1).tolist()

    res_repl = compute_results(labels, truth)
    res_repl["time"] = [elapsed_density + elapsed_phm]
    res_repl["time_density"] = [elapsed_density]
    res_repl["time_phm"] = [elapsed_phm]
    res_repl["n_components"] = [len(clusters)]
    try:
        save_json(fname, res_repl)
    except Exception as exc:
        ## Clean up
        Path(fname).unlink()
        raise exc

    return res_repl
