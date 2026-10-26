import numpy as np
from pathlib import Path
from time import time
from pathos.multiprocessing import ProcessPool

from helpers import init_results, load_json, save_json, compute_results
from config import leiden_resolution, hdbscan_mcs, REPLICATIONS, CORES, part1_base
from model_fitting import fit_predict_mclust2, fit_predict_sc3s, fit_predict_leiden, fit_predict_kmeans, fit_predict_hdbscan #, fit_predict_skgmm


results_kmeans = init_results()
results_gmmRand = init_results()
results_gmmKmeans = init_results()
results_gmmKpp = init_results()
results_gmmHier = init_results()
results_sc3s5 = init_results()
results_sc3s50 = init_results()
results_leiden = init_results()
results_hdbscan = init_results()

def loop_func_baseline_methods(rs, fname_base, fp_func, args):
    fname = fname_base.format(rpl=rs)
    if Path(fname).exists():
        return load_json(fname)

    start = time()
    labels = fp_func(X, random_state=rs, **args)
    elapsed = time() - start
    
    if isinstance(labels, np.ndarray):
        labels = labels.tolist()

    res_repl = compute_results(labels)
    res_repl["time"] = [elapsed]

    try:
        save_json(fname, res_repl)
    except Exception as exc:
        ## Clean up
        Path(fname).unlink()

    return res_repl

targets_vec = [
        # (part1_base.format(mthd="gmmRand", rpl="{rpl}"), fit_predict_skgmm, results_gmmRand, {"n_clusters": 10, "init": "random_from_data"}),
        # (part1_base.format(mthd="gmmKmeans", rpl="{rpl}"), fit_predict_skgmm, results_gmmKmeans, {"n_clusters": 10, "init": "kmeans"}),
        # (part1_base.format(mthd="gmmK++", rpl="{rpl}"), fit_predict_skgmm, results_gmmKpp, {"n_clusters": 10, "init": "k-means++"}),
        # (part1_base.format(mthd="sc3s5", rpl="{rpl}"), fit_predict_sc3s, results_sc3s5, {"n_clusters": 10, "n_runs": 5}),
        (part1_base.format(mthd="mclust", rpl="{rpl}"), fit_predict_mclust2, results_gmmHier, {"n_clusters": 10}),
        (part1_base.format(mthd="sc3s50", rpl="{rpl}"), fit_predict_sc3s, results_sc3s50, {"n_clusters": 10, "n_runs": 50}),
        (part1_base.format(mthd="leiden", rpl="{rpl}"), fit_predict_leiden, results_leiden, {"resolution": leiden_resolution}),
        (part1_base.format(mthd="kmeans", rpl="{rpl}"), fit_predict_kmeans, results_kmeans, {"n_clusters": 10}),
        (part1_base.format(mthd="hdbscan", rpl="{rpl}"), fit_predict_hdbscan, results_hdbscan, {"mcs": hdbscan_mcs})
    ]

random_states = [x+1 for x in range(REPLICATIONS)]


for tgt in targets_vec:
    print("Fitting for {}".format(tgt[0]))
    fname_base = [tgt[0]] * REPLICATIONS
    fp_func = [tgt[1]] * REPLICATIONS
    args = [tgt[3]] * REPLICATIONS
    res_obj = tgt[2]
    idx = 0
    with ProcessPool(nodes=CORES) as pool:
        for res in pool.uimap(loop_func_baseline_methods, random_states, fname_base, fp_func, args):
            idx += 1
            print("\tReplicate: {}".format(idx), end="\r")
            for key in res_obj.keys():
                res_obj[key].extend(res.get(key, []))

        print("\n")

print("Complete!")