import numpy as np
from pathlib import Path
from time import time

from helpers import init_results, load_json, save_json, compute_results
from config import leiden_resolution, hdbscan_mcs, REPLICATIONS, CORES, part1_base
from model_fitting import fit_predict_mclust2, fit_predict_sc3s, fit_predict_leiden, fit_predict_kmeans, fit_predict_hdbscan #, fit_predict_skgmm
from model_fitting import ProcessPool

results_kmeans = init_results()
results_mclust = init_results()
results_sc3s50 = init_results()
results_leiden = init_results()
results_hdbscan = init_results()

def loop_func_baseline_methods(X, truth, rs, fname_base, fp_func, args):
    fname = fname_base.format(rpl=rs)
    if Path(fname).exists():
        return load_json(fname)

    start = time()
    labels = fp_func(X, random_state=rs, **args)
    elapsed = time() - start
    
    if isinstance(labels, np.ndarray):
        labels = labels.tolist()

    res_repl = compute_results(labels, truth)
    res_repl["time"] = [elapsed]

    try:
        save_json(fname, res_repl)
    except Exception as exc:
        ## Clean up
        Path(fname).unlink()

    return res_repl

targets_vec = [
        (part1_base.format(mthd="mclust", rpl="{rpl}"), fit_predict_mclust2, results_mclust, {"n_clusters": 10}),
        (part1_base.format(mthd="sc3s50", rpl="{rpl}"), fit_predict_sc3s, results_sc3s50, {"n_clusters": 10, "n_runs": 50}),
        (part1_base.format(mthd="leiden", rpl="{rpl}"), fit_predict_leiden, results_leiden, {"resolution": leiden_resolution}),
        (part1_base.format(mthd="kmeans", rpl="{rpl}"), fit_predict_kmeans, results_kmeans, {"n_clusters": 10}),
        (part1_base.format(mthd="hdbscan", rpl="{rpl}"), fit_predict_hdbscan, results_hdbscan, {"mcs": hdbscan_mcs})
    ]
