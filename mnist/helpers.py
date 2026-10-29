import json
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.cluster import contingency_matrix

#######################
## Wrappers for JSON ##
#######################
def save_json(filename, obj):
    with open(filename, 'w') as handle:
        json.dump(obj, handle, indent=4)

def load_json(filename):
    with open(filename, 'r') as handle:
        b = json.load(handle)
        
    return b

##################################
## For printing in multiprocess ##
##################################
def lock_print(lock, *args, **kwargs):
    with lock:
        print(*args, **kwargs)

    return

#################################
## Results Tracking/Evaluation ##
#################################
def init_results():
    return {
        "NMI": [], "ARI": [], "Acc": [], "predicted": [], 
        "models": [], 
        "time": [], "time_density": [], "time_phm": [],
        "n_components": []
    }


def compute_results(pred, actual):
    res = init_results()
    res["NMI"].append(normalized_mutual_info_score(actual, pred))
    res["ARI"].append(adjusted_rand_score(actual, pred))
    res["Acc"].append(cluster_accuracy(actual, pred))
    res["predicted"].append(pred)
    return res

def results_print(results_list, digits=3):
    res = np.round([np.mean(results_list), np.median(results_list), np.min(results_list), np.max(results_list)], digits)
    print(
        "Mean: {}".format(res[0]),
        "Median: {}".format(res[1]),
        "Range: ({}, {})".format(res[2], res[3]),
        sep="\n", end="\n\n"
    )

def gen_results_dict(results, mthd, q=[0.05, 0.95]):
    out = {
        "Method": mthd.replace("\n", " ")
    }
    # out["ARI"] = np.round(np.median(results["ARI"]), 3)
    out["ARI"] = np.round(np.mean(results["ARI"]), 3)
    if len(results["ARI"]) > 1:
        # out["ARI Range"] = "({:.2f}, {:.2f})".format( *np.quantile(results["ARI"], q))
        out["ARI SD"] = np.round(np.std(results["ARI"]), 3)
        if out["ARI SD"] == 0:
            out["ARI SD"] = "---"
    else:
        out["ARI SD"] = "---"

    # out["NMI"] = np.round(np.median(results["NMI"]), 3)
    out["NMI"] = np.round(np.mean(results["NMI"]), 3)
    if len(results["NMI"]) > 1:
        # out["NMI Range"] = "({:.2f}, {:.2f})".format( *np.quantile(results["NMI"], q, method="inverted_cdf") )
        out["NMI SD"] = np.round(np.std(results["NMI"]), 3)
        if out["NMI SD"] == 0:
            out["NMI SD"] = "---"
    else:
        out["NMI SD"] = "---"


    # out["Acc"] = np.round(np.median(results["Acc"]), 3)
    out["Acc"] = np.round(np.mean(results["Acc"]), 3)
    if len(results["Acc"]) > 1:
        # out["Acc Range"] = "({:.2f}, {:.2f})".format( *np.quantile(acc, q, method="inverted_cdf") )
        out["Acc SD"] = np.round(np.std(results["Acc"]), 3)
        if out["Acc SD"] == 0:
            out["Acc SD"] = "---"
    else:
        out["Acc SD"] = "---"
    
    return out

def cluster_accuracy(y_true, y_pred):
    cm = contingency_matrix(y_true, y_pred)
    row_ind, col_ind = linear_sum_assignment(-cm)
    return cm[row_ind, col_ind].sum() / len(y_true)