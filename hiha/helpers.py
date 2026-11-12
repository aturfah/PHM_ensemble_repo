import numpy as np
import plotnine as p9
import pandas as pd
from functools import reduce

from PHMCluster import PHMCluster
from distbn_helpers import mvgmm_generator

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

def load_params(filename):
    readRDS_func = robjects.r['readRDS']
    params = readRDS_func(filename)
    return convert_R_params_object2(params)

def convert_R_params_object2(res, probscale=1):
    clusters = []
    densities = []
    with localconverter(robjects.default_converter + numpy2ri.converter):
        for clust_params in res:
            ## Extract params into dictionary numpy arrays
            param_dict = {}
            for val in clust_params.items():
                param_dict[val.name] = val.value

            ## Construct PHMClusters for each mixture component
            G = param_dict["prob"].shape[0]
            for comp_idx in range(G):
                prob = param_dict["prob"][comp_idx]
                mean = np.reshape(param_dict["mean"][:, comp_idx], (param_dict["mean"].shape[0], ))
                var = np.reshape(param_dict["var"][:, :, comp_idx], (param_dict["var"].shape[0], param_dict["var"].shape[1]))
                density, sampler = mvgmm_generator([prob * probscale], [mean], [var])
                clust = PHMCluster([prob * probscale], [sampler], [density])
                clusters.append(clust)

    return clusters, densities

def fuzzy_contingency(class_labels, posterior_matrix, scale_rows=True):
    unique_classes, class_indices = np.unique(
        class_labels, return_inverse=True
    )
    N, K = posterior_matrix.shape
    C = len(unique_classes)

    table = np.empty((K, C), dtype=posterior_matrix.dtype)

    for k in range(K):
        table[k, :] = np.bincount(
            class_indices, weights=posterior_matrix[:, k], minlength=C
        )
    if scale_rows:
        row_sums = table.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        table = table / row_sums

    return pd.DataFrame(table, columns=unique_classes)
