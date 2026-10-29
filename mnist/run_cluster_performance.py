import numpy as np
import pandas as pd
import plotnine as p9
from time import time
from multiprocessing import Manager
from functools import partial

from config import REPLICATIONS, CORES, images_directory
from helpers import init_results, gen_results_dict
from model_fitting import ProcessPool

from baseline import targets_vec, loop_func_baseline_methods, results_kmeans, results_mclust, results_sc3s50, results_leiden, results_hdbscan
from phm_single import loop_func_phm_single
from phm_ensemble_global import loop_func_ensemble_global
from phm_ensemble_localized import loop_func_ensemble_localized

from prepare_data import load_tsne, load_labels

proj_tsne = load_tsne()
sampleY_t = load_labels()

########################
### Baseline Methods ###
########################
print("==== Running Baseline Methods ====", end="\r")
random_states = [x+1 for x in range(REPLICATIONS)]
for tgt in targets_vec:
    print("\nFitting for {}".format(tgt[0]))
    fname_base = [tgt[0]] * REPLICATIONS
    fp_func = [tgt[1]] * REPLICATIONS
    args = [tgt[3]] * REPLICATIONS
    res_obj = tgt[2]
    idx = 0

    loop_func = partial(loop_func_baseline_methods, proj_tsne, sampleY_t)
    
    with ProcessPool(nodes=CORES) as pool:
        for res in pool.uimap(loop_func, random_states, fname_base, fp_func, args):
            idx += 1
            print("\tReplicate: {}".format(idx), end="\r")
            for key in res_obj.keys():
                res_obj[key].extend(res.get(key, []))

print("\n==== Completed Baseline ====\n")


##############################
### PHM (Single Replicate) ###
##############################
print("\n==== Running Single-Replicate PHM ====")
results_phm_single = init_results()
random_states = [x for x in range(1, 101)]
pfx_list = REPLICATIONS * ["mclustALL"]
idx = 0
with ProcessPool(nodes=CORES) as pool:
    loop_func = partial(loop_func_phm_single, proj_tsne, sampleY_t)
    for res in pool.uimap(loop_func, random_states, pfx_list):
        idx += 1
        print("Completed: {}".format(idx), end="\r")
        if res is None:
            continue

        for key in results_phm_single.keys():
            results_phm_single[key].extend(res.get(key, []))

print("\n==== Finished Single-Replicate PHM ====\n")


#############################
### PHM-Ensemble (Global) ###
#############################
print("\n==== Running Ensemble PHM (Global) ====")
pfx = "mclustSubaggALL{num_pcs}_r{idx}"
random_states = [x for x in range(1, 100 + 1)]
results_phm_ens_global = init_results()
counter = 0
m = Manager()
print_lock = m.Lock()

with ProcessPool(nodes=CORES) as pool:
    lock_list = [print_lock] * REPLICATIONS
    pfx_list = [pfx] * REPLICATIONS
    loop_func = partial(loop_func_ensemble_global, proj_tsne, sampleY_t)
    for res in pool.uimap(loop_func, random_states, lock_list, pfx_list):
        counter += 1
        for key in results_phm_ens_global.keys():
            results_phm_ens_global[key].extend(res.get(key, []))

        print("Completed: {}".format(counter), flush=True, end="\r")

print("\n==== Finished Ensemble PHM (Global) ====\n")


################################
### PHM-Ensemble (Localized) ###
################################
print("\n==== Running Ensemble PHM (Localized) ====")
random_states = [x for x in range(REPLICATIONS)]
counter = 0
start = time()
results_phm_ens_local = init_results()
with ProcessPool(nodes=CORES) as pool:
    loop_func = partial(loop_func_ensemble_localized, proj_tsne, sampleY_t)
    for res in pool.uimap(loop_func, random_states):
        counter += 1
        print("Completed: {} | {:.2f} minutes".format(counter, (time() - start)/60 ), flush=True, end="\r")
        for key in results_phm_ens_local.keys():
            results_phm_ens_local[key].extend(res.get(key, []))

print("\n==== Finished Ensemble PHM (Localized) ====\n")



################################
### Baseline Methods Results ###
################################
text_kmeans = "k-means"
text_leiden = "Leiden"
text_sc3s = "SC3s"
text_hdbscan = "HDBSCAN"
text_hclust = "Hier. Clust"
text_mclust = "mclust"
text_phm_single = "PHM"
text_phm_ens_global = "PHM-Ensemble\n(Global)"
text_phm_ens_local = "PHM-Ensemble\n(Localized)"
text_oracle = "Oracle"

## Oracle and Hierarchical Clustering come from mnist_density_estimation_v2.R
results_dict_hclust = {"Method": text_hclust, "ARI": 0.789, "ARI SD": "---", "NMI": 0.840, "NMI SD": "---", "Acc": 0.888, "Acc SD": "---"}
results_dict_oracle = {"Method": text_oracle, "ARI": 0.938, "ARI SD": "---", "NMI": 0.925, "NMI SD": "---", "Acc": 0.971, "Acc SD": "---"}

## Results in Table Format 
gen_results_dict_lambda = lambda x, y: gen_results_dict(x, y, [0, 1])
results_df = pd.DataFrame([
    gen_results_dict_lambda(results_kmeans, text_kmeans),
    gen_results_dict_lambda(results_mclust, text_mclust),
    gen_results_dict_lambda(results_sc3s50, text_sc3s),
    gen_results_dict_lambda(results_leiden, text_leiden),
    gen_results_dict_lambda(results_hdbscan, text_hdbscan),
    gen_results_dict_lambda(results_phm_single, text_phm_single),
    gen_results_dict_lambda(results_phm_ens_global, text_phm_ens_global.replace("\n", " ")),
    gen_results_dict_lambda(results_phm_ens_local, text_phm_ens_local.replace("\n", " ")),
    results_dict_hclust,
    results_dict_oracle
])


print("\n==== Method Performance Table ====\n")
results_df = results_df.sort_values("ARI")
print(results_df.to_markdown(index=False), end="\n\n")

## Boxplot for visualizing the range
data_dict = {
    text_phm_ens_global: results_phm_ens_global["ARI"],
    text_phm_ens_local: results_phm_ens_local["ARI"],
    text_phm_single: results_phm_single["ARI"],
    text_hclust: [results_dict_hclust["ARI"]],
    text_leiden: results_leiden["ARI"],
    text_mclust: results_mclust["ARI"],
    text_hdbscan: results_hdbscan["ARI"],
    text_sc3s: results_sc3s50["ARI"],
    text_kmeans: results_kmeans["ARI"]
}

boxplot_df = pd.DataFrame([
    {"Method": name, "Value": val} for name, vals in data_dict.items()
    for val in vals
])

boxplot_df["Method"] = pd.Categorical(
    boxplot_df["Method"], 
    categories=reversed([
        text_leiden, 
        text_kmeans, 
        text_hclust,
        text_sc3s, 
        text_mclust,
        text_hdbscan,
        text_phm_single,
        text_phm_ens_global,
        text_phm_ens_local]), 
    ordered=True)

plot_baseline = (
    p9.ggplot(boxplot_df, p9.aes(x='Method', y='Value', fill='Method')) +
    p9.scale_y_continuous(breaks=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]) +
    p9.geom_boxplot() +
    p9.theme_bw() +
    p9.geom_hline(yintercept=results_dict_oracle["ARI"], linetype="dashed") +
    p9.scale_fill_brewer(type="qual", palette="Paired") +
    p9.coord_flip() +
    p9.labs( 
        x='Method', 
        y='ARI'
    ) +
    p9.theme(figure_size=(10, 5), text=p9.element_text(size=10))
)

plot_baseline.save(images_directory.joinpath("mnist_baseline_results.png"))