from config import REPLICATIONS, CORES
from helpers import init_results
from model_fitting import ProcessPool
from time import time
from multiprocessing import Manager
from functools import partial

from baseline import targets_vec, loop_func_baseline_methods
from phm_single import loop_func_phm_single
from phm_ensemble_global import loop_func_ensemble_global
from phm_ensemble_localized import loop_func_ensemble_localized

from prepare_data import proj_tsne, sampleY_t

########################
### Baseline Methods ###
########################
random_states = [x+1 for x in range(REPLICATIONS)]

for tgt in targets_vec:
    print("Fitting for {}".format(tgt[0]))
    fname_base = [tgt[0]] * REPLICATIONS
    fp_func = [tgt[1]] * REPLICATIONS
    args = [tgt[3]] * REPLICATIONS
    res_obj = tgt[2]
    idx = 0
    loop_func_partial = partial(loop_func_baseline_methods, X=proj_tsne)

    with ProcessPool(nodes=CORES) as pool:
        for res in pool.uimap(loop_func_baseline_methods, random_states, fname_base, fp_func, args):
            idx += 1
            print("\tReplicate: {}".format(idx), end="\r")
            for key in res_obj.keys():
                res_obj[key].extend(res.get(key, []))

        print("\n")

print("==== Completed Baseline ====")


raise RuntimeError("AAAAHHH")


##############################
### PHM (Single Replicate) ###
##############################
results_phmGMMHier = init_results()
random_states = [x for x in range(1, 101)]
pfx_list = REPLICATIONS * ["mclustALL"]
idx = 0
with ProcessPool(nodes=CORES) as pool:
    for res in pool.uimap(loop_func_phm_single, random_states, pfx_list):
        idx += 1
        print("Completed: {}".format(idx), end="\r")
        if res is None:
            continue

        for key in results_phmGMMHier.keys():
            results_phmGMMHier[key].extend(res.get(key, []))

print("\n==== Finished Single-Replicate PHM ====")


#############################
### PHM-Ensemble (Global) ###
#############################

pfx = "mclustSubaggALL{num_pcs}_r{idx}"
random_states = [x for x in range(1, 100 + 1)]
results_phmAuto = init_results()
counter = 0
m = Manager()
print_lock = m.Lock()

with ProcessPool(nodes=CORES) as pool:
    lock_list = [print_lock] * REPLICATIONS
    pfx_list = [pfx] * REPLICATIONS
    for res in pool.uimap(loop_func_ensemble_global, random_states, lock_list, pfx_list):
        counter += 1
        for key in results_phmAuto.keys():
            results_phmAuto[key].extend(res.get(key, []))

        print("Completed: {}".format(counter), flush=True, end="\r")

print("\n==== Finished Ensemble PHM (Global) ====")
################################
### PHM-Ensemble (Localized) ###
################################

random_states = [x for x in range(REPLICATIONS)]

counter = 0
start = time()
results_phmSubaggLeiden = init_results()
with ProcessPool(nodes=CORES) as pool:
    for res in pool.uimap(loop_func_ensemble_localized, random_states):
        counter += 1
        print("Completed: {} | {:.2f} minutes".format(counter, (time() - start)/60 ), flush=True, end="\r")
        for key in results_phmSubaggLeiden.keys():
            results_phmSubaggLeiden[key].extend(res.get(key, []))

print("\n==== Finished Ensemble PHM (Localized) ====")

